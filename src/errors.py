import numpy as np

def RMS(test, groundtruth, nb_dim=3):
    """
    Compute the Root Mean Square error between two 1D arrays.
    nb_dim = 0 pour X, 1 pour Y, 2 pour Z
    """
    assert len(test) == len(groundtruth), "Les deux tableaux doivent avoir la même taille"
    if (nb_dim == 0):
        Xtest = test.T[0]
        XGT = groundtruth.T[0]
        return np.sqrt(np.sum((Xtest-XGT)**2)/len(Xtest))
    elif (nb_dim == 1):
        Ytest = test.T[1]
        YGT = groundtruth.T[1]
        return np.sqrt(np.sum((Ytest-YGT)**2)/len(Ytest))
    elif (nb_dim == 2):
        Ztest = test.T[2]
        ZGT = groundtruth.T[2]
        return np.sqrt(np.sum((Ztest-ZGT)**2)/len(Ztest))
    else:
        return np.sqrt(np.sum((test-groundtruth)**2)/len(test)) #TODO:  à verifier

    
    