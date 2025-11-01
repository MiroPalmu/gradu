# Beam-plasma derivation

Starting from relativistic Vlasov equation (in the Guassian units):

\begin{align}
\partial_t f_j
+ \mathbf{v} \cdot \nabla f_j
+ q_j ( \mathbf{E} + \boldsymbol{\beta} \times \mathbf{B}) \cdot \nabla_\mathbf{p} f_j = 0
\end{align}

our goal is to derive continuity equation:

\begin{align}
\partial_t n_j + \nabla \cdot (n_j \mathbf{v}_j) = 0
\end{align}

and the Euler equation for cold plasma:

\begin{align}
\partial_t \mathbf{p}_j
+ (\mathbf{v}_j \cdot \nabla) \mathbf{p}_j
= q_j (\mathbf{E} + \boldsymbol{\beta}_j \times \mathbf{B})
\end{align}

where by cold plasma we mean that
f$_j = n_j(\mathbf{x}, t)\delta(\mathbf{p} - \mathbf{p}_j(\mathbf{x}, t))$.
