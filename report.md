# Compte rendu des réunions
12/02

Organisation générale:
- réunion hebdomadaire de 1h-1h30, à priori le jeudi
    - à moi d'ouvrir la réunion: ce que j'ai fait / compte faire / mes pistes
    - trace écrite sur git
- possibilité (plus tard) de 2-3h de travail en commun avec Yoann

Idée générale du projet:
Ajouter des infos aux graphes de Bruijn (dBG) (colorés ou pas ?) pour permettre de reconstruire les haplotypes d'origine. Le but n'est pas d'avoir une reconstruction rapide mais un graphe compact en espace mémoire (donc d'ajouter un minimum d'infos, et les infos les plus pertinentes possibles).

À faire (par complexité ~ croissante):
- résumé lecture d'articles: ce que j'en retire
- pistes à explorer: ajout d'infos au dBG, mètriques pour estimer à quel point on reconstruit les haplotypes d'origines
- installation d'outils existants (pggb, Minigraph Cactus, ggcat) et test sur des exemples simples (génome de moustiuque/ un seul chromosome)
- prise en main nextflow

__________________
fait durant la semaine 10/02-19/02:
- formation sécurité
- lecture d'articles: intérêt et utilisation des pangénomes; comparaison de méthodes; implémentation de méthodes (pggb, gccat, mac-dBG)
- lecture Rust book + tutos
- installation et utilisation de certains outils:
    - données = moustiques tigres Aalf3 et Aalf5 (Aalf5: assemblage de 3 chromosomes + 1493 scafolds non placés; Aalf3: 574 scafolds avec indication de chromosome)
- utilisation de wfmash pour partitionner les séquences par chromosomes: mappage de Aalf3 (query) sur Aalf5 (target)
    - on retrouve bien 3 clusters principaux (6 en réalités mais 3 ne contiennent que 2 séquences chacun)
    - grosso-modo, les scaffolds de Aalf3 sont mappés sur le bon chromosomes de Aalf5, mais pas tous -> partitions selon wfmash ou selon les indications des fastas ? dans le premier cas, on fait quoi des scaffolds de Aalf5 non placés ? dans l'autre, des scaffolds de Aalf3 mal placés ?
- utilisation de pggb sur chaque communauté: pertinance de créer un graphe de variation sur un génome assemblé en scaffolds -> mapping des scaffolds ?

questions:
- pour Aalf3: numérotation chr[1-3].XX -> XX arbitraire ou signification ? d'ailleurs, comment est obtenu le chromosome ? par analyse post-séquençage et mapping ?