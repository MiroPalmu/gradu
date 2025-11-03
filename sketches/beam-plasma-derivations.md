# Beam-plasma derivation

Starting from relativistic Vlasov equation (in the Guassian units):

\begin{align}
\partial_t f_\alpha
+ \mathbf{v} \cdot \nabla f_\alpha
+ q_\alpha ( \mathbf{E} + \boldsymbol{\beta} \times \mathbf{B}) \cdot \nabla_\mathbf{p} f_\alpha = 0
\end{align}

our goal is to derive continuity equation:

\begin{align}
\partial_t n_\alpha + \nabla \cdot (n_\alpha \mathbf{v}_\alpha) = 0
\end{align}

and the Euler equation for cold plasma:

\begin{align}
\partial_t \mathbf{p}_\alpha
+ (\mathbf{v}_\alpha \cdot \nabla) \mathbf{p}_\alpha
= q_\alpha (\mathbf{E} + \boldsymbol{\beta}_\alpha \times \mathbf{B})
\end{align}

where by cold plasma we mean that
$f_\alpha = n_\alpha(\mathbf{x}, t)\delta(\mathbf{p} - \mathbf{p}_\alpha(\mathbf{x}, t))$.


## Moments

We want to evaluate the moments $M(\mathbf{p})$
(we drop subscripts $\alpha$ for clarity):

\begin{align}
\int M \left(
\partial_t f
+ \mathbf{v} \cdot \nabla f
+ q ( \mathbf{E} + \boldsymbol{\beta} \times \mathbf{B}) \cdot \nabla_\mathbf{p} f = 0
\right) d\mathbf{p}
\end{align}


## First term

\begin{align}
\int M \partial_t f d\mathbf{p} = \partial_t \int M f d\mathbf{p}
\end{align}


## Second term

\begin{align}
\int M \mathbf{v} \cdot \nabla f d\mathbf{p}
= \partial_i \int M v^i f d\mathbf{p}
\end{align}


## Third term

First define $a^i = qE^i + \frac{q}{c}\epsilon^i_{lk}v^l B^k$
and calculate:

\begin{align}
\partial_{p_i} a^i = \frac{q}{c}\epsilon^i_{lk}\partial_{p_i}v^l B^k
\end{align}

using:

\begin{align}
\partial_{p_i} v^l
&= \partial_{p_i} \frac{p^l}{m \gamma} \\
&= \frac{\delta^l_i}{m \gamma} -  \frac{p^l}{m\gamma^2}\partial_{p_i} \gamma \\
&= \frac{\delta^l_i}{m \gamma} -  \frac{p^l}{m\gamma^2}\partial_{p_i} \sqrt{1 + \frac{p^np_n}{m^2c^2}} \\
&= \frac{\delta^l_i}{m \gamma} -  \frac{p^lp_i}{m^3c^2\gamma^3}
\end{align}

\begin{align}
\therefore \partial_{p_i} a^i
&= \frac{q}{c}\epsilon^i_{lk}\left(\frac{\delta^l_i}{m \gamma} -  \frac{p^lp_i}{m^3c^2\gamma^3}\right) B^k \\
&= 0
\end{align}

\begin{align}
\int M \mathbf{a} \cdot \nabla_\mathbf{p} f d\mathbf{p}
&= \sum_i \int \int M a^i \partial_{p_i} f dp_i \prod_{i \neq j}dp_j \\
&= - \sum_i \int \int \partial_{p_i}( M a^i ) f dp_i \prod_{i \neq j}dp_j \\
&= - \sum_i \int \int (a^i \partial_{p_i} M + M \partial_{p_i} a^i) f dp_i \prod_{i \neq j}dp_j \\
&= - \int a^i \partial_{p_i} M f d\mathbf{p} \\
\end{align}


# Overall

\begin{align}
\partial_t \int M f_\alpha d\mathbf{p}
+ \partial_i \int M v^i f_\alpha d\mathbf{p}
= \int a^i \partial_{p_i} M f_\alpha d\mathbf{p}
\end{align}

Now for $f_\alpha = n_\alpha(\mathbf{x}, t)\delta(\mathbf{p} - \mathbf{p}_\alpha(\mathbf{x}, t))$
we get continuity equation with $M = 1$:

\begin{align}
\partial_t n_\alpha + \partial_i ( v^i_\alpha n_\alpha ) = 0
\end{align}

and Euler equation for cold plasma with $M = p^k$ after some manipulation:

\begin{align}
\partial_t (p^k_\alpha n_\alpha)
+ \partial_i (p^k_\alpha v^i_\alpha n_\alpha)
&= \int a^i \partial_{p_i} p^k f_\alpha d\mathbf{p} \\
&= a^k(\mathbf{p}_\alpha) n_\alpha \\
&= q_\alpha (E^k + \frac{\epsilon^k_{mn}v^m_\alpha B^n}{c}) n_\alpha
\end{align}

Insert continuity equation to the LHS:

\begin{align}
\partial_t (p^k_\alpha n_\alpha) + \partial_i (p^k_\alpha v^i_\alpha n_\alpha)
&= \partial_t p^k_\alpha n_\alpha + p^k_\alpha \partial_t n_\alpha
+ \partial_i p^k_\alpha v^i_\alpha n_\alpha
+ p^k_\alpha \partial_i (v^i_\alpha n_\alpha) \\
&= \partial_t p^k_\alpha n_\alpha - p^k_\alpha \partial_i(v^i_\alpha n_\alpha)
+ \partial_i p^k_\alpha v^i_\alpha n_\alpha
+ p^k_\alpha \partial_i (v^i_\alpha n_\alpha) \\
&= \partial_t p^k_\alpha n_\alpha + \partial_i p^k_\alpha v^i_\alpha n_\alpha
\end{align}

So overall:

\begin{align}
\partial_t \mathbf{p}_\alpha n_\alpha + n_\alpha (\mathbf{v}_\alpha \cdot \nabla) \mathbf{p}_\alpha
&= q_\alpha (\mathbf{E} + \boldsymbol{\beta}_\alpha \times \mathbf{B}) n_\alpha \\
\partial_t \mathbf{p}_\alpha + (\mathbf{v}_\alpha \cdot \nabla) \mathbf{p}_\alpha
&= q_\alpha (\mathbf{E} + \boldsymbol{\beta}_\alpha \times \mathbf{B})
\end{align}