# Freezing detection from thermographic videos (Winke et al. 2022, XXX)

This is a Matlab script used to detect freezing episodes from thremographic camera videos (Testo885) on hot-plate experiments (used in Winke et al. 2022, XXX).


- Define ROI that includes the hotplate area.
- Analyze video.

Briefly, the algorithm:
- Analize videos on the blue channel (better image substraction with mice in our range of background temperatures).
- Binarize individual image frames based on automatic thresholding Otsu's method (Ref).
- Mask the area corresponded to mouse for individual image frames.
- Perform consecutive frame image substraction and compute instantaneous % area change.
- 