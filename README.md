# Electrostatics and magnetostatics on the hypersphere

While studying at Masaryk University in Brno, I wrote my bachelor thesis - Electrostatics and magnetostatics on the hypersphere. The title says it all. My task was to solve scalar and vector Poisson equations for electrostatic and magnetostatic potentials. I managed to find a solution for the electrostatic potential. Still, unfortunately, the magnetostatic (vector) potential was too much of a task and would require extra time, which I did not have at the time. I managed to defend my thesis with the grade A in June of 2024. If you are interested in the Thesis itself you can have a look here [An Internal Link](/electrostatics_and_magnetostatics_on_the_hypersphere_valentik.pdf)







Materials and code that I use for my bachelor thesis - electrostatics and magnetostatics on a hypersphere

Scripts that are in a folder called "scripts" ( :) ) are unfortunately not well documented (sometimes not at all).

Most of my scripts are used to decompose functions on the hypersphere (sphere in four dimensions) into eigenfunctions of the Laplace operator. In other words hyperspherical harmonics (the most interesting are the Gegenbauer polynomials that are used).

I also have scripts that allow me to visualize the functions on the hypersphere. My script allows me to plot any hyperspherical function to 3d space using hyperstereographic projection. The resulting image can look something like this:

<img src="https://github.com/Valentyk/Thesis/assets/146948734/04bc89ee-d1e3-4fa3-8c78-9ae125e5dffb" width="300" align="right">

The function shown is the function of $\zeta$ where $\zeta$ is one of the hyperspehrical coordinates ($\phi$, $\theta$, $\zeta$).

The plot in the picture carries some kind of visible numerical artifacts, especially in the yellow parts of the plot. If I have enough time for my thesis I will try to get rid of them.
