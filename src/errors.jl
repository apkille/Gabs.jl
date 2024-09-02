
STATE_ERROR = lazy"Either the input covariance matrix is not
    square or its dimensions do not match the mean vector. An
    arbitrary N-mode Gaussian state is characterized by a mean
    vector of length 2N and a covariance matrix of size 2N x 2N."

CHANNEL_ERROR = lazy"Either an input matrix is not square or 
    there is a dimension mismatch. An arbitrary N-mode Gaussian 
    channel is characterized by a displacement vector of length
    2N, a transform matrix of size 2N x 2N, and a noise matrix
    of size 2N x 2N." 

ACTION_ERROR = lazy"The number of modes for the Gaussian channel
    does not match the number of modes for the Gaussian state."
