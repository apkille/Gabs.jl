
STATE_ERROR = lazy"Either the input covariance matrix is not
    square or its dimensions do not match the mean vector. An
    arbitrary N-mode Gaussian state is characterized by a mean
    vector of length 2N and a covariance matrix of size 2N x 2N."

UNITARY_ERROR = lazy"Either the input symplectic matrix is not square 
    or there is a dimension mismatch. An arbitrary N-mode Gaussian 
    unitary operator is characterized by a displacement vector of 
    length 2N and a symplectic matrix of size 2N x 2N." 

CHANNEL_ERROR = lazy"Either an input matrix is not square or 
    there is a dimension mismatch. An arbitrary N-mode Gaussian 
    channel is characterized by a displacement vector of length
    2N, a transform matrix of size 2N x 2N, and a noise matrix
    of size 2N x 2N." 

ACTION_ERROR = lazy"The number of modes for the Gaussian operator
    does not match the number of modes for the Gaussian state."

GENERALDYNE_ERROR = lazy"The number of modes for the measurement state
    does not match the number of modes for the projected subsytem."

WIGNER_ERROR = lazy"The length of your input array does not match the
    number of modes for the Gaussian state."

HEATMAP_ERROR = lazy"The input Gaussian state describes more than one mode. 
    A heat map visualization for a multi-mode Gaussian state is not possible."
