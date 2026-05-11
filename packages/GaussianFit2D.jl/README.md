# STM Molecular Pattern Fit

Reconnaissance exploratoire de motifs moléculaires dans des images STM Nanonis
`.sxm`.

Ce projet remplace l'idée trop fragile “sweep global de N gaussiennes + BIC” par
un pipeline plus adapté aux images STM :

1. lecture du fichier `.sxm` ;
2. conversion des unités (`m → nm`, `A → pA`) ;
3. prétraitement STM : soustraction de plan et correction ligne-par-ligne ;
4. fusion multi-vues alignées (`Z fwd/bwd` par défaut) ;
5. détection généreuse de blobs/features candidates ;
6. extraction de chaînes moléculaires géométriquement cohérentes ;
7. raffinement local par gaussiennes 2D contraintes ;
8. export des features/chaînes moléculaires et plot de diagnostic.

Le fit gaussien n'est donc plus utilisé comme vérité globale, mais comme
raffinement des features détectées.

Deux modes objectifs de sélection de `N` sont disponibles :

- `--robust-sweep` : diagnostic avancé par mélange libre de `N` gaussiennes 2D
  dans la ROI. Ce mode ne doit pas être interprété comme le nombre moléculaire :
  une composante libre peut modéliser une aile, une asymétrie ou l'enveloppe.
- `--chain-sweep` : chaîne 2D ordonnée, plus comparable au programme 1D : les
  centres sont contraints le long de l'axe principal de la molécule, avec
  espacements, décalages latéraux et largeurs bornés. Par défaut le fit se fait
  dans un tube autour de l'axe, pas sur toute la ROI.

## Utilisation

```bash
cd /home/durif/Git/GaussianFit2D.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Afficher les canaux disponibles :

```bash
julia --project=. run.jl path/to/your_image.sxm --info
```

Lancer une reconnaissance avec plot :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --channel Z --direction fwd \
  --stride 2 --contrast auto --fusion-channels Z \
  --threshold-sigma 2.5 --min-distance-px 10 \
  --chain-min-length 3 --chain-max-spacing-nm 1.6 --chain-max-angle-deg 65 \
  --output-dir results/251206_034_pattern
```

Balayer le seuil de bruit pour vérifier la stabilité de la chaîne :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --threshold-sweep 2.0:0.5:4.0 --fusion-channels Z \
  --chain-max-spacing-nm 1.6 --chain-max-angle-deg 65 \
  --output-dir results/251206_034_pattern
```

Inspecter toutes les images contenues dans le `.sxm` :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --inspect-images --output-dir results/251206_034_all_images
```

Sélection robuste par mélange libre de gaussiennes 2D :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --robust-sweep --robust-n-min 0 --robust-n-max 10 \
  --multistart 8 --cv-folds 4 --stride 3 \
  --output-dir results/251206_034_robust_sweep
```

Sélection par chaîne 2D ordonnée :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --chain-sweep --chain-n-min 2 --chain-n-max 10 \
  --chain-spacing-min-nm 0.35 --chain-spacing-max-nm 0.75 \
  --chain-lateral-max-nm 0.35 --chain-fit-width-nm 0.45 \
  --multistart 8 --cv-folds 4 --stride 3 \
  --output-dir results/251206_034_chain_sweep
```

Traiter tous les couples canal/direction (`Z fwd`, `Z bwd`, `Current fwd`, etc.) :

```bash
julia --project=. run.jl path/to/your_image.sxm \
  --all-images --max-features 50 --max-fit-features 8 \
  --output-dir results/251206_034_all_images
```

Sorties :

- `molecular_features.tsv` : amplitude, position et largeur de chaque feature ;
- `molecular_chains.tsv` : chaînes acceptées, score, espacement moyen, régularité ;
- `molecular_pattern_overview.png` : image brute, image prétraitée + détections,
  modèle raffiné, résidus.
- `all_images_overview.png` / `all_images_summary.tsv` en mode multi-images.
- `robust_model_selection.tsv`, `robust_selected_lobes.tsv`,
  `robust_model_selection.png` en mode `--robust-sweep`.
- `chain_model_selection.tsv`, `chain_selected_lobes.tsv`, `chain_axis.tsv`,
  `chain_axis_profile.tsv`, `chain_model_selection.png` en mode
  `--chain-sweep`.

Les coordonnées exportées sont des coordonnées relatives à l'image
(`0 → SCAN_RANGE`) en nm, pas des coordonnées absolues de scène Nanonis.
Les scans backward sont retournés horizontalement afin d'être alignés
spatialement avec les scans forward.

## Options importantes

- `--flatten plane+rows` : recommandé pour STM. Alternatives : `none`, `plane`,
  `rows`.
- `--contrast auto|bright|dark` : choisir le signe attendu du motif.
- `--threshold-sigma` : seuil de détection en unités de bruit robuste.
- `--min-distance-px` : distance minimale entre features candidates.
- `--max-features` : nombre maximal de features moléculaires à raffiner.
- `--max-fit-features` : limite les candidats effectivement passés au fit
  non-linéaire, pour éviter des fits trop lents/instables sur les images Current.
- `--no-fit` : seulement détection, sans raffinement non-linéaire.
- `--fusion-channels Z|Current|Z,Current|all` : canaux utilisés pour la carte
  d'évidence commune. Par défaut `Z`, c'est-à-dire topographie fwd+bwd, car le
  courant peut contenir plus de structures parasites.
- `--chain-min-length`, `--chain-min/max-spacing-nm`, `--chain-max-angle-deg`,
  `--min-chain-score` : contraintes structurelles qui transforment des blobs en
  chaînes moléculaires acceptées.
- `--max-path-branches` : largeur de recherche du meilleur chemin dans le graphe
  de candidats. Augmenter si la chaîne semble tronquée, diminuer si les chemins
  deviennent instables.
- `--chain-sweep` : fit une chaîne de gaussiennes 2D ordonnées. Les centres sont
  paramétrés par `t` le long de l'axe principal et `u` perpendiculairement :
  `centre_k = origine + t_k axe + u_k axe_perp`.
- `--chain-spacing-min/max-nm` : bornes sur `t_{k+1}-t_k`. Pour comparer au fit
  1D, utiliser des bornes proches de celles du modèle 1D.
- `--chain-lateral-max-nm` : tolérance latérale maximale autour de l'axe.
- `--chain-fit-width-nm` : demi-largeur du tube de pixels utilisé par le fit. À
  diminuer pour se rapprocher d'une slice 1D ; à augmenter pour fitter plus de
  l'enveloppe 2D.
- `--chain-support-threshold-fraction`, `--chain-support-padding-nm` : détectent
  automatiquement la partie longitudinale active du profil avant le fit.
- `--chain-sigma-parallel/perp-min/max-nm` : bornes des largeurs gaussiennes
  parallèles/perpendiculaires à la chaîne. Rappel : `FWHM = 2.3548 * sigma`.

## Interprétation de `N`

`--robust-sweep` et `--chain-sweep` ne répondent pas exactement à la même
question :

- le mélange libre compte les blobs 2D statistiquement utiles dans la ROI ;
- la chaîne contrainte compte les positions le long d'un axe moléculaire, même si
  certains maxima sont faibles ou partiellement fusionnés.

Si le profil 1D donne `N=6`, le test comparable côté 2D est donc : même fichier
SXM, même axe, même région longitudinale, tube étroit, mêmes bornes
d'espacement/FWHM. Le `robust-sweep N=6` libre n'est pas cette comparaison.

## Limites assumées

La prochaine étape pour un résultat vraiment physique est d'ajouter un **template
moléculaire** : positions relatives des lobes/features, rotation, translation et
éventuelles petites déformations. Le pipeline actuel prépare cette étape en
fournissant des features robustes, plutôt qu'un BIC global non fiable.
