# Compte rendu des réunions

### Réunion du 12/02

Organisation générale:
- réunion hebdomadaire de 1h-1h30, à priori le jeudi
    - à moi d'ouvrir la réunion: ce que j'ai fait / compte faire / mes pistes
    - trace écrite sur git
- possibilité (plus tard) de 2-3h de travail en commun avec Yoann

Idée générale du projet:
Ajouter des infos aux graphes de Bruijn (dBG) (colorés ou pas ?) pour permettre de reconstruire les haplotypes d'origine. Le but n'est pas d'avoir une reconstruction rapide mais un graphe compact en espace mémoire (donc d'ajouter un minimum d'infos, et les infos les plus pertinentes possibles).

**À faire (par complexité ~ croissante):**
- résumé lecture d'articles: ce que j'en retire
- pistes à explorer: ajout d'infos au dBG, mètriques pour estimer à quel point on reconstruit les haplotypes d'origines
- installation d'outils existants (pggb, Minigraph Cactus, ggcat) et test sur des exemples simples (génome de moustiuque/ un seul chromosome)
- prise en main nextflow

__________________
### Fait durant la semaine 10/02-19/02:
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
### Réunion du 20/02

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
    - regarder le fonctionnement de smoothxg
    - séance whiteboard avec Francesco + exploration perso

- moyen/moyen-long terme:
    - implémenter une métrique, la tester avec 1 puis 2 haplotypes
    - ajouter des infos pour améliorer le score