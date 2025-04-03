# Compte rendu des réunions

## Réunion 12/02

Organisation générale:
- réunion hebdomadaire de 1h-1h30, à priori le jeudi
    - à moi d'ouvrir la réunion: ce que j'ai fait / compte faire / mes pistes
    - trace écrite sur git
- possibilité (plus tard) de 2-3h de travail en commun avec Yoann

Idée générale du projet:
Ajouter des infos aux graphes de Bruijn (dBG) (colorés ou pas ?) pour permettre de reconstruire les haplotypes d'origine. Le but n'est pas d'avoir une reconstruction rapide mais un graphe compact en espace mémoire (donc d'ajouter un minimum d'infos, et les infos les plus pertinentes possibles).

**À faire (par complexité ~ croissante):**
- [X] résumé lecture d'articles: ce que j'en retire
- [-] pistes à explorer: ajout d'infos au dBG, mètriques pour estimer à quel point on reconstruit les haplotypes d'origines
- [X] installation d'outils existants (pggb, Minigraph Cactus, ggcat) et test sur des exemples simples (génome de moustiuque/ un seul chromosome)
- [ ] prise en main nextflow

__________________
### 10/02-19/02:
- formation sécurité
- lecture d'articles: intérêt et utilisation des pangénomes; comparaison de méthodes; implémentation de méthodes (pggb, gccat, mac-dBG)
- lecture Rust book + tutos
- installation et test de pggb, ggcat,...
- manipulation données: moustiques tigres Aalf3 et Aalf5 (Aalf5: assemblage de 3 chromosomes + 1493 scafolds non placés; Aalf3: 574 scafolds avec indication de chromosome)
- utilisation de wfmash pour partitionner les séquences par chromosomes: mappage de Aalf3 (query) sur Aalf5 (target)
    - on retrouve bien 3 clusters principaux (6 en réalités mais 3 ne contiennent que 2 séquences chacun)
    - grosso-modo, les scaffolds de Aalf3 sont mappés sur le bon chromosomes de Aalf5, mais on a certaines erreurs
    - pertinence ? 97% des scaffolds sont déjà placés -> mieux vaut utiliser les indications du fasta

__________________
## Réunion 20/02

Réflexion sur les métriques associés au dBG:
- nb de noeud
- degrés des noeuds
- comptage de kmers
- nb de nouveaux noeuds créés en ajoutants un haplotype
- edit distance: nb de modifications (SNP, indel,...) pour passer d'une branche à une autre

Plusieurs échelles possibles
- graphe
- haplotype
- fenêtre, de taille variable
- noeud, et voisinage

Propriétés importantes:
- définie rigoureusement -> résultat parfaitement reproductible
- calculable -> applicable en pratique
- maintenabilité dynamique (dynamic graph algorithms): mise à jour possible après ajout ou retrait de noeuds/arrêtes -> permet des optimisations locales (cf descente de gradient)

Autre:
- algo de compression des données, et des graphes en particuliers (gzip, )

**À faire:**  
- court terme:
    - [X] regarder le fonctionnement de smoothxg
    - [X] séance whiteboard avec Francesco + exploration perso

- moyen/moyen-long terme:
    - [ ] implémenter une métrique, la tester avec 1 puis 2 haplotypes
    - [ ] ajouter des infos pour améliorer le score

__________________
### 21/02 au 27/02
- lecture documentation:
    - smoothxg
    - odgi
    - graphAligner: fonctionne avec des graphes biorientés avec ou sans overlapp (VG, dBG, assembly graph,...) -> "alignment graph"

- format des graphes:
    - Bifrost -> dBG sous forme de gfa, C++ API for graph manipulation
    - ggcat -> juste le set de kmers (peut encore être comprimé avec matchtig et eulertig)

__________________
## Réunion 27/02
- comparer la taille des fichiers permettant / ne permettant pas de reconstruire un haplotype. Dans les deux cas, prendre la version la plus épurée, mais non compressée possible
- pour mesurer la reconstructibilité des haplotypes: regarder le nb de "breakpoints" le long d'un haplotype
    - pour les gfa: 0
    - pour les dBG: beaucoup (↘ avec k, ↗ avec le nb de genomes)

__________________
### 28/02 au 06/03
- début code rust: recherche crate pour manipuler les kmers; structure basique
- comparaison taille fichier (texte) unitigs/matchtigs/vg

__________________
### 07/03 au 21/03
- article sshash: hash pour la query efficiente des kmers. pas d'implémentation rust -> on reste sur du plus basique
- poursuite code rust: prototype qui construit le graphe et parcourt les haplotypes

__________________
## Réunion 21/03
- mêtriques à regarder:
    - % de breackpoint, fenêtres glissantes
    - % de edges non utilisés
- 3 origines (a minima) aux breackpoints:
    - SNP et compagnie -> comparer le nb de edges en ajoutant un haplo
    - répétitions
    - (k-1) overlap non présent dans les séquence -> comparer graphes pour k et k+1
- lever les indéterminations: quel infos ajouter ?
    - checkpoints -> Djikistra
    - railway switches
    - forbidden intersections  
L'intérêt des checkpoints: méthode mixte, hiérachique. on peut utiliser différentes techniques pour différentes portions

__________________
### 21/03 au 01/04
- début implementation Breadth-first search:
    - note: un même kmer (canonique) peut apparaitre deux fois dans le plus court chemin (une fois dans chaque sens)
    - comment trouver le point d'arrivée optimal ? première idée: recherche dichotomique -> beaucoup trop long
    - deuxième idée: lors du parcours en largeur, on énumère en même temps les kmers de l'haplotype, et on vérifie qu'ils sont bein dans la bonne profondeur
    - troisième idée: notion de direction ? -> permettrait d'utiliser l'algo A*
__________________
## Réunion 01/04
- mêtrique: **nb de breackpoints supplémentaires en ajoutant d'autres haplotypes**: exemple avec plusieurs haplo d'un chromosome humain
- plus court chemin:
    - Djikistra par les deux bouts -> réduit la complexité temporelle
    - recherche dichotomique pour allonger un chemin au maximum

__________________
### 02/04 au ???
suite BFS: ça marche pas. pour split de 1 à 16, chemin trouvé seulement pour 1, 5 et 13

## Comparaison crates rust
#### Ragnar Groot: packed-seq, ptr-hash, simd-minimizers
+crates très optimisées  
-super pointu -> dur d'utilisation  
-bcp de crates différentes -> pb de compatibilités entre elles; certaines seulement sur git et pas crates.io

#### debruijn
+du consortium 10X genomics; beaucoup d'utilisateurs  
+pas mal de fonctionnalités rassemblés en un seul lieu
-pas possible d'itérer directement sur les kmer canoniques

#### rust-bio
+very well documented and widely used  
+-lots of functionalities, in a wide range of applications, but nothing very specific

#### needletail: for fasta parsing and kmer streaming
+easy to use  
+slightly faster than bio reader