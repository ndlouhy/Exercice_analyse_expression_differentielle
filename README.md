# Exercice_analyse_expression_differentielle

## Introduction
Ce dépôt contient le travail réalisé dans le cadre de l’exercice de sélection proposé par Simon Ducheix, chercheur au sein de l'unité de recherche de l'institut du thorax à Nantes et Yuna Blum chercheuse de l'Institut de Génétique & Développement de Rennes.
L’objectif de cet exercice est d’explorer un jeu de données de transcriptomique (comptages de gènes) réalisé sur 26 échantillons appartenant à trois groupes biologiques distincts :

Groupe 1 : échantillons 1-8\
Groupe 2 : échantillons 9-17\
Groupe 3 : échantillons 18-26\

## Objectif
L'objectif de l’analyse est double. Dans un premier temps, il faut déterminer la structure des données et évaluer si elles sont exploitable. Dans un second temps il faut réaliser une analyse d'expression différentielle afin d'identifier les tops gènes les plus différentiellement exprimés entre les conditions.

## Analyse

### Import du jeu de données
Dans un premier temps, il faut importer le jeu de données.
<pre> 
file_path <- "data/gene_count.xls"

data_file <- read.table(file_path, header = T, sep = "\t")
</pre>

Description :
Le jeu de données présente 35 colonnes pour 33 808 lignes.
Il est présenté ainsi :
Col 1 : gene_id : identifiant Ensembl du gène.
Col 2-27 : Sample_1-26 : Colonnes de comptages de lectures pour chaque échantillons.
Col 28 : gene_name : Nom du gène
Col 29 : gene_chr : Chromosome ou se trouve le gène.
Col 30 : gene_start : Position de départ du gène.
Col 31 : gene_end : Position de fin du gène.
Col 32 : gene_strand : Brin du chromosome ou si trouve le gène.
Col 33 : gene_length : Longueur du gène.
Col 34 : gene_biotype : Catégorie fonctionnelle du gène
Col 35 : gene_description : Description de la fonction du gène.
Col 36 : tf_family : Catégorie de famille de facteurs de transcription du gène.
										
<img width="3745" height="25" alt="image" src="https://github.com/user-attachments/assets/49272cb1-7e14-4848-b249-bb75fb8253f7" />

## Résultats

## Conclusion
