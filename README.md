# Electrostatics and magnetostatics on the hypersphere

While studying at Masaryk University in Brno, I wrote my bachelor thesis - Electrostatics and magnetostatics on the hypersphere. The title says it all. My task was to solve scalar and vector Poisson equations for electrostatic and magnetostatic potentials. I managed to find a solution for the electrostatic potential. Still, unfortunately, the magnetostatic (vector) potential was too much of a task and would require extra time, which I did not have at the time. I managed to defend my thesis with the grade A in June of 2024. If you are interested in the Thesis, you can download it [here](https://is.muni.cz/th/z48op/?lang=en) or [here](/electrostatics_and_magnetostatics_on_the_hypersphere_valentik.pdf).

I wrote some code for my bachelor thesis which I share here for anyone interested.

## Hypersphere

What do I mean by hypersphere? In my thesis by hypersphere I mean 3-sphere $\mathbb{S}^3$ - three dimensional sphere imbedded in $\mathbb{R}^4$. We can define hypersphere $\mathbb{S}^3$ like this

$$\mathbb{S}^3 = \\{ (x_1,x_2,x_3,x_4) \in \mathbb{R}^4 \vert \sqrt{x_1^2 + x_2^2 + x_3^2 + x_4^2} = r \\},$$

where $(x_1,x_2,x_3,x_4)$ is a point on $\mathbb{R}^4$ (expressed by Cartesian coordinates) and $r$ is the radius of the hypersphere. To make things little bit easier for me I set the hypersphere radius $r = 1.$ 

### Hyperspherical coordinates

Cartesian coordinates are not the best coordinates for describing the hypersphere (trust me...) and more "natrual" choice of coordinates can make the work with the hypersphere a lot easier. One of the more "natural" coordinate system that I used in my thesis is hyperspherical coordinate system represented by the map $\alpha$

$$
\alpha: \qquad \begin{array}{l}
            x_1 = \cos\varphi, \\
            x_2 = \sin\varphi\cos\vartheta, \\
            x_3 = \sin\varphi\sin\vartheta\cos\zeta, \\
            x_4 = \sin\varphi\sin\vartheta\sin\zeta.
        \end{array}
$$

where $\varphi \in (0,\pi)$, $\vartheta \in (0,\pi)$ and $\zeta \in (0,2\pi)$ are the new coordinates - angles - of the hypersphere ($\zeta$ is analogous to azimuthal angle in spherical coordinates). We can derive the inverse map $\alpha^{-1}$

$$
\alpha^{-1}: \qquad \begin{array}{l}
        \varphi = \arctan\frac{\sqrt{x_2^2+x_3^2+x_4^2}}{x_1}+\frac{\pi}{2}, \\
        \vartheta = 
            \begin{cases}
            \arctan\frac{\sqrt{x_3^2+x_4^2}}{x_2}+\frac{\pi}{2}, & \mbox{if } x_2 \neq 0, x_3 \neq 0 \mbox{ and } x_4 \neq 0, \\
            \text{undefined}, & \mbox{for } x_2 = x_3 = x_4 = 0,
            \end{cases} \\
        \zeta = 
            \begin{cases} 
                \arctan\left(\frac{x_4}{x_3}\right), & \mbox{if } x_3>0, \\ 
                \arctan\left(\frac{x_4}{x_3}\right) + \pi, & \mbox{if } x_3<0 \mbox{ and } x_4 \ge 0, \\
                \arctan\left(\frac{x_4}{x_3}\right) - \pi, & \mbox{if } x_3<0 \mbox{ and } x_4<0, \\
                +\frac{\pi}{2}, & \mbox{if } x_3=0 \mbox{ and } x_4>0, \\
                -\frac{\pi}{2}, & \mbox{if } x_3=0 \mbox{ and } x_4<0, \\
                \text{undefined}, & \mbox{if } x_3 = x_4=0.
            \end{cases}
    \end{array}
$$

which will be useful later :).

### Stereographic projection

For the human mind, the hypersphere $\mathbb{S}^3$ is not easily imaginable. Fortunately we can map the hypersphere $\mathbb{S}^3$ onto the $\mathbb{R}^3$ by stereographic map (more on stereographic projection [here](https://en.wikipedia.org/wiki/Stereographic_projection)) $\beta$

$$
\beta: \qquad \begin{array}{l}
        \mu = \frac{x_1}{1-x_4}, \\
        \nu = \frac{x_2}{1-x_4}, \\
        \xi = \frac{x_3}{1-x_4},
    \end{array}
$$

where $\mathbf{p} = (x_1, x_2, x_3, x_4)$ represents a point on the hypersphere (and therefore it must be true that $\sqrt{x_1^2 + x_2^2 + x_3^2 + x_4^2} = 1$) and $\mu$, $\nu$ and $\xi$ are the new coordinates of $\mathbf{p}$ on $\mathbb{R}^3$. Inverse stereographic projection $\beta^{-1}$ will prove to be useful so let me write it down here

$$
\beta^{-1}: \qquad \begin{array}{l}
        x' = \frac{2\mu}{\mu^2+\nu^2+\xi^2+1}, \\
        y' = \frac{2\nu}{\mu^2+\nu^2+\xi^2+1}, \\
        z' = \frac{2\xi}{\mu^2+\nu^2+\xi^2+1}, \\
        w' = \frac{\mu^2+\nu^2+\xi^2-1}{\mu^2+\nu^2+\xi^2+1}.
        \end{array}
$$

### Visualisation of functions on the hypersphere - finally some code!

Now we can visualize the functions on the hypersphere $f(\varphi, \vartheta, \zeta)$ using the stereographic projection. How? Let us start with the grid of $N^3$ (we will call $N$ the 'grid number') points in the closed subspace $M = \\{p = (x_1, x_2, x_3) \vert x_1 \in [-2,2], x_2 \in [-2,2], x_3 \in [-2,2]\\}$ of $\mathbb{R}^3$. In the code, this process looks like this

```python 
N_fig = 120                                              # Grid number defines the "resolution" of the Figures

x_values = np.linspace(-2,2,N_fig)                       # Defiing the grid for the plotting
y_values = np.linspace(-2,2,N_fig)
z_values = np.linspace(-2,2,N_fig)
x,y,z = np.meshgrid(x_values, y_values, z_values)                               
```

Then, using the inverse map $\beta^{-1}$ (inverse stereographic projection), we will express each point of the grid in $R^4$. Finally, we can use the inverse hyperspherical map $\alpha^{-1}$, we will end up with the hyperspherical coordinates of each point of the grid (from our subspace $M$) into which we want to project the function $f(\varphi, \vartheta, \zeta)$. These inverse maps are defined in the code by the functions `R3_R4` and `R4_S3` like this

``` python
def R3_R4(x_values, y_values, z_values):                 # Inverse stereographic projection
    X = (2*x_values)/(x_values**2+y_values**2+z_values**2+1) 
    Y = (2*y_values)/(x_values**2+y_values**2+z_values**2+1) 
    Z = (2*z_values)/(x_values**2+y_values**2+z_values**2+1) 
    W = (x_values**2+y_values**2+z_values**2-1)/(x_values**2+y_values**2+z_values**2+1) 
    return(X,Y,Z,W)

def R4_S3(X,Y,Z,W):                                      # Inverse hyperspherical map
    zeta = np.arctan2(W,Z)
    theta = np.arctan2(np.sqrt(Z**2+W**2),Y)
    phi = np.arctan2(np.sqrt(W**2+Z**2+Y**2),X)
    return(phi, theta, zeta)
```

To confuse the reader - the code denotes the coordinates in $\mathbb{R}^4$ as `X`, `Y`, `Z`, `W` and in $\mathbb{R}^3$ as `x_values`, `y_values`, `z_values` instead of $x_1$, $x_2$, $x_3$, $x_4$ or $\mu$, $\nu$, $\xi$ respectively.

Now, we can assign the hyperspherical coordinate to each point of the grid in the subspace $M$ of $\mathbb{R}^3$ and that means that we are also able to assign the function value of the function $f(\varphi, \vartheta, \zeta)$. The hard part is over... uff...

Now the only thing that remains to do is to plot the function. I decided to use the isosurfaces (surfaces with constant function value) to be plotted and used [Plotly](https://plotly.com/python/) graphing library to plot the figure of the isosurfaces of the function.
```python
func_values_g = np.real(func_g(zeta_g,theta_g,phi_g))    # Evaluation of the func_g for each point of the plot grid 

    fig = go.Figure(layout=layout, data = go.Volume(     # Making the plot
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value = func_values_g.flatten(),
        isomin = np.min(func_values_g),
        isomax = np.max(func_values_g),
        # isomin = -1,
        # isomax = 1,
        surface_count = 10,
        opacity = 0.5,
        #surface_fill = 0.9,
        colorscale = "plasma",
        caps = dict(x_show = False, y_show = False, z_show = False)
    ))

    fig.update_layout(scene_camera=camera)

    fig.show()
```

In my thesis I chose the function $f = \cos\varphi \cos\vartheta \cos\zeta$ as an example of the visualization and the output image of my script looks like this
<img src="https://github.com/user-attachments/assets/d151408b-4107-434d-a6af-36673b11841e" width="300" align="right">


Let us have the function $f(\varphi, \vartheta, \zeta)$ (notice that the function's arguments are the hyperspherical angles). With the use of the map $\alpha^{-1}$ transform the function $f(\varphi, \vartheta, \zeta)$ with the   



## Junk, don't mind this part :)

Materials and code that I use for my bachelor thesis - electrostatics and magnetostatics on a hypersphere

Scripts that are in a folder called "scripts" ( :) ) are unfortunately not well documented (sometimes not at all).

Most of my scripts are used to decompose functions on the hypersphere (sphere in four dimensions) into eigenfunctions of the Laplace operator. In other words hyperspherical harmonics (the most interesting are the Gegenbauer polynomials that are used).

I also have scripts that allow me to visualize the functions on the hypersphere. My script allows me to plot any hyperspherical function to 3d space using hyperstereographic projection. The resulting image can look something like this:

<img src="https://github.com/Valentyk/Thesis/assets/146948734/04bc89ee-d1e3-4fa3-8c78-9ae125e5dffb" width="300" align="right">

The function shown is the function of $\zeta$ where $\zeta$ is one of the hyperspehrical coordinates ($\phi$, $\theta$, $\zeta$).

The plot in the picture carries some kind of visible numerical artifacts, especially in the yellow parts of the plot. If I have enough time for my thesis I will try to get rid of them.
