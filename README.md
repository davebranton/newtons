# Dependencies

You'll need `libpng-dev`

# Build

Build with

```
g++ -std=c++11 -o NewstonsFractal -O3 main.cpp -lm -lpng -lstdc++ -lpthread
```

# Run

It's got a few parameters.

```
Usage: NewtonsFractal [-w value] [-h value] -o value [-d value] [-r value] [-i value] [-s value] [-p value] [-f value] [-n value] [-m value] [-t value] [-v value] [-e value] [-a value] [-l value] [-u value] [-y value]
        -w/--width               : Image Width
        -h/--height              : Image Height
        -o/--file                : Output filename
        -d/--degree              : Degree of polynomial
        -r/--real_ellipse_factor : Flatten to ellipse along real axis
        -i/--imag_ellipse_factor : Flatten to ellipse along imaginary axis
        -s/--scale               : Scale image
        -p/--spiral_scale        : Spiral scale start
        -f/--spiral_factor       : Spiral scale factor
        -n/--n_threads           : Number of threads to use
        -m/--max_iters           : Max iterations
        -t/--rotation            : Radians to rotate roots by
        -v/--overlap             : Root circle overlap factor
        -e/--edges               : Colour only the edges
        -a/--saturation          : Colour saturation
        -l/--value               : Colour value
        -u/--start_hue           : Start Hue
        -y/--hue_multiplier      : Hue multiplier
```