---
title: "Individual-based disease models"
excerpt: "Simulations to be parameterized to mesocosm data to further explore diversity-disease relationship"
layout: single
classes: wide
entries_layout: grid
header:
  teaser: images/mesocosm/spread_map_still.jpg
gallery:
  - url: images/mesocosm/spread_map_detsub2.gif
    image_path: images/mesocosm/spread_map_detsub2.gif
    alt: "simulated disease spread in 2 species communities"
  - url: images/mesocosm/spread_map_detsub6.gif
    image_path: images/mesocosm/spread_map_detsub6.gif
    alt: "simulated disease spread in 2 species communities"
---

In efforts to futher explore how community composition affects disease risk, I wrote an individual-based disease model to accompany my mesocosm experiment. My data from the mesocosm experiment will allow me to validate my model so I can explore additional starting parameters. The model is a spatially explicit discrete time model that tracks state changes for all individuals over the duration of the epidemic. Once I figure out how to properly render the mathematical notation on this website, I'll post it. Until then, check out simulation results for 2-species and 6-species communities:  

{% include gallery caption="Simulated disease spread in 2- and 6-species communities. Colors represent different host species and black dots indicated infected individuals. Model assumes species vary in competency (ability to acquire and transmit pathogen), transmission decays with distance, and transmission rates vary with time." %}