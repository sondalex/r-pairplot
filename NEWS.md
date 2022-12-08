# pairplot 0.1.2

## Changes

* Changed `pairplot::pair_geom_histogram` scaling based on neighbour internals.
  Before the bounds for min max scaling of the data were set to
  the neighbour `panel_scales_x` and `panel_scales_y` attributes.
  In some case this would lead
  to having negative lower bounds, in particular
  for `panel_scales_y`. y-axis scales often have a negative offset.
  In those cases the lower bound of the scales is negative.
  Now, the bounds are set to the transformed data of neighbour plots.
  **It might still be unstable**. I would have to look more at internals
  of `ggplot2` library.


# pairplot 0.1.1

## Changes

* `diag_share_ylim` is deprecated in favor of `diag_share_lim`.
* Remove pairgrid unused arguments.
* Add possible knowledge of surrounding subplots for diagonal subplots. Works also with corner plot.

## Bug Fixes

* **Important**: Fixed min-max scaling error in `pair_geom_histogram`.


# pairplot 0.1.0

* The layout system relies on [patchwork](https://patchwork.data-imaginist.com/index.html)
* Now accessible on GitHub.
