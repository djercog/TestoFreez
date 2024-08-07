# Freezing detection from thermographic videos.

This repository contains Matlab scripts used to detect freezing episodes from thremographic camera videos (Testo885) on hot-plate experiments (used in Winke, Aby, Jercog et al. 2024, XXX).


- Define ROI that includes the hotplate area.
- Analyze video.

Briefly, the video analysis algorithm:
- Extract blue channel from RGB videos (better image substraction with mice in our range of background temperatures).
- Binarize individual image frames based on automatic thresholding Otsu's method (Ref).
- Mask the "mouse area" for individual image frames.
- Perform consecutive frame image substraction and compute instantaneous % mouse area change.
- Threshold % area change to detect freezing periods.

![image](https://user-images.githubusercontent.com/28762337/186142723-75874576-17e5-428c-9239-4fd8a2f7d2af.png)
