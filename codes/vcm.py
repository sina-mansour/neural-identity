# Python implementation of morphometricity (variance component modelling)
# Author of python implementation: Sina Mansour L.
# contact: sina.mansour.lakouraj@gmail.com

# implementing the codes from sabuncu et. al.
# original source code is available in Matlab
# http://people.csail.mit.edu/msabuncu/morphometricity/Morphometricity.m

import numpy as np
import scipy.linalg as la


# Matlab translations
def _mrdivide(B, A):  # return B/A
    return la.lstsq(A.T, B.T)[0].T


def _mldivide(B, A):  # return B\A
    return la.solve(B, A)


def _variance_component_model_reimplemented(y, X, K, alg=0, tol=1e-4, MaxIter=100, verbal=False):
    '''
    % Input -
    % y is an Nsubj x 1 vector of phenotypes (trait values)
    % X is an Nsubj x Ncov covariate matrix (that contains "nuisance variables
    % such as age, sex, site dummy, etc)
    % K is an Nsubj x Nsubj anatomical similarity matrix (ASM)
    %  K has to be a symmetric, positive semi-definite matrix with its diagonal
    %  elements averaging to 1.
    %   If K is not non-negative definite, its negative eigenvalues will be set to zero,
    %   and a warning will be printed to the command window
    % alg is the algorithm for the ReML; default alg = 0
    %   alg = 0 - use the average information
    %   alg = 1 - use the expected information (Fisher scoring)
    %   alg = 2 - use the observed information
    % tol is the tolerance for the convergence of the ReML algorithm; default tol = 1e-4
    % MaxIter is the maximum number of iterations for the ReML algorithm; default MaxIter = 100
    %
    % Output - a dictionary containing the following information:
    % flag indicates the convergence of the ReML algorithm
    %   flag = 1 - the ReML algorithm has converged
    %   flag = 0 - the ReML algorithm did not converged
    % m2 is the morphometricity estimate
    % SE is the standard error of the morphometricity estimate
    % Va is the total anatomical/morphological variance
    % Ve is the residual variance
    % Lnew is the ReML likelihood when the algorithm is converged
    '''

    # NOTE: in python D is a vector but in matlab it is a diagonal matrix
    # [U, D] = eig(K);
    D, U = la.eigh(K)

    # check whether the GRM is non-negative definite
    if (np.min(D) < 0):
        print('WARNING: the GRM is not non-negative definite! Set negative eigenvalues to zero')
        D[D < 0] = 0   # set negative eigenvalues to zero
        K = _mrdivide(np.dot(U, np.diag(D)), U)   # reconstruct the GRM (K = U*D/U;)

    # derived quantities
    Nsubj = len(y)   # calculate the total number of subjects
    Vp = np.var(y)   # calculate the phenotypic variance

    # --- initialization

    # initialize the anatomical variance and residual variance
    Va = Vp / 2
    Ve = Vp / 2

    # initialize the covariance matrix
    V = np.dot(Va, K) + np.dot(Ve, np.eye(Nsubj))

    # initialize the projection matrix (P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;)
    P = _mrdivide(
        (
            np.eye(Nsubj) -
            np.dot(
                _mrdivide(
                    (_mldivide(V, X)),
                    np.dot(
                        _mrdivide(X.T, V),
                        X
                    )
                ),
                X.T
            )
        ),
        V
    )

    # --- EM algorithm

    if verbal:
        print('---------- EM algorithm ----------')

    # use the expectation maximization (EM) algorithm as an initial update

    # update the anatomical variance
    # Va = (Va^2*y'*P*K*P*y+trace(Va*eye(Nsubj)-Va^2*P*K))/Nsubj;
    Va = float((
        np.dot(
            np.dot(
                np.dot(
                    np.dot(
                        ((Va**2) * y.T),
                        P
                    ),
                    K
                ),
                P
            ),
            y
        ) +
        np.trace(
            (Va * np.eye(Nsubj)) -
            np.dot(
                ((Va**2) * P),
                K
            )
        )) / Nsubj
    )

    # update the residual variance
    # Ve = (Ve^2*y'*P*P*y+trace(Ve*eye(Nsubj)-Ve^2*P))/Nsubj;   % update the residual variance
    Ve = float((
        np.dot(
            np.dot(
                np.dot(
                    ((Ve**2) * y.T),
                    P
                ),
                P
            ),
            y
        ) +
        np.trace(
            np.dot(
                Ve,
                np.eye(Nsubj)
            ) -
            ((Ve**2) * P)
        )) / Nsubj
    )

    # set negative estimates of the variance component parameters to Vp*1e-6
    if (Va < 0):
        Va = (10**(-6)) * Vp
    if (Ve < 0):
        Ve = (10**(-6)) * Vp

    # update the covariance matrix
    # V = Va*K+Ve*eye(Nsubj);
    V = (Va * K) + (Ve * np.eye(Nsubj))

    # update the projection matrix
    # P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V
    P = _mrdivide(
        (
            np.eye(Nsubj) -
            np.dot(
                _mrdivide(
                    _mldivide(V, X),
                    np.dot(
                        _mrdivide(X.T, V),
                        X
                    )
                ),
                X.T
            )
        ),
        V
    )

    # calculate the log determinant of the covariance matrix
    E = la.eigvals(V)
    logdetV = np.sum(np.log(np.real(E)))

    # initialize the ReML likelihood
    Lold = np.inf
    # Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y
    Lnew = (
        (-0.5 * logdetV) -
        (0.5 * np.log(la.det(np.dot(_mrdivide(X.T, V), X)))) -
        (0.5 * np.dot(np.dot(y.T, P), y))
    )

    # --- ReML

    if verbal:
        print('---------- ReML iterations ----------')

    # initialize the total number of iterations
    iter = 0

    # criteria of termination
    while ((abs(Lnew - Lold) >= tol) and (iter < MaxIter)):

        # new iteration
        iter += 1
        Lold = Lnew

        if verbal:
            print('---------- REML Iteration- {} ----------'.format(iter))

        # update the first-order derivative of the ReML likelihood

        # score equation of the anatomical variance
        # Sg = -1/2*trace(P*K)+1/2*y'*P*K*P*y
        Sg = float((-0.5 * np.trace(np.dot(P, K))) + (0.5 * np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), y)))

        # score equation of the residual variance
        # Se = -1/2*trace(P)+1/2*y'*P*P*y
        Se = ((-0.5 * np.trace(P)) + (0.5 * np.dot(np.dot(np.dot(y.T, P), P), y)))

        # construct the score vector
        S = np.matrix([Sg, Se]).T

        # update the information matrix
        if alg == 0:
            # average information
            # Igg = 1/2*y'*P*K*P*K*P*y;
            Igg = float(0.5 * np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), K), P), y))

            # Ige = 1/2*y'*P*K*P*P*y;
            Ige = float(0.5 * np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), P), y))

            # Iee = 1/2*y'*P*P*P*y;
            Iee = (0.5 * np.dot(np.dot(np.dot(np.dot(y.T, P), P), P), y))

        elif alg == 1:
            # expected information
            # Igg = 1/2*trace(P*K*P*K)
            Igg = float(0.5 * np.trace(np.dot(np.dot(np.dot(P, K), P), K)))

            # Ige = 1/2*trace(P*K*P)
            Ige = float(0.5 * np.trace(np.dot(np.dot(P, K), P)))

            # Iee = 1/2*trace(P*P)
            Iee = (0.5 * np.trace(np.dot(P, P)))

        elif alg == 2:
            # observed information
            # Igg = -1/2*trace(P*K*P*K)+y'*P*K*P*K*P*y
            Igg = float(0.5 * np.trace(np.dot(np.dot(np.dot(P, K), P), K))) + np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), K), P), y)

            # Ige = -1/2*trace(P*K*P)+y'*P*K*P*P*y
            Ige = float(0.5 * np.trace(np.dot(np.dot(P, K), P))) + np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), P), y)

            # Iee = -1/2*trace(P*P)+y'*P*P*P*y;
            Iee = -(0.5 * np.trace(np.dot(P, P))) + np.dot(np.dot(np.dot(np.dot(y.T, P), P), P), y)

        # construct the information matrix
        I = np.matrix([[Igg, Ige], [Ige, Iee]])

        # update the variance component parameters
        T = np.matrix([Va, Ve]).T + _mldivide(I, S)
        Va = T[0, 0]
        Ve = T[1, 0]

        # set negative estimates of the variance component parameters to Vp*1e-6
        if (Va < 0):
            Va = (10**(-6)) * Vp
        if (Ve < 0):
            Ve = (10**(-6)) * Vp

        # update the covariance matrix
        V = (Va * K) + (Ve * np.eye(Nsubj))

        # update the projection matrix
        # P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V
        P = _mrdivide(
            (
                np.eye(Nsubj) -
                np.dot(
                    _mrdivide(
                        _mldivide(V, X),
                        np.dot(
                            _mrdivide(X.T, V),
                            X
                        )
                    ),
                    X.T
                )
            ),
            V
        )

        if (np.isnan(V).any()):
        # if (np.isnan(V).any() or np.isinf(Lnew)):
            flag = 0
            m2 = np.nan
            SE = np.nan
            Va = np.nan
            Ve = np.nan
            Lnew = np.nan
            return {'flag': flag, 'm2': m2, 'SE': SE, 'Va': Va, 'Ve': Ve, 'Lnew': Lnew, 'dist': np.nan}

        # calculate the log determinant of the covariance matrix
        E = la.eigvals(V)
        logdetV = np.sum(np.log(np.real(E + np.finfo(E.dtype).eps)))

        # update the ReML likelihood
        # Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y
        Lnew = (
            (-0.5 * logdetV) -
            (0.5 * np.log(la.det(np.dot(_mrdivide(X.T, V), X)))) -
            (0.5 * np.dot(np.dot(y.T, P), y))
        )

    # --- morphometricity estimate and standard error

    # morphometricity estimate
    m2 = Va / (Va + Ve)

    # update the information matrix at the final estimates
    if alg == 0:
        # average information
        # Igg = 1/2*y'*P*K*P*K*P*y;
        Igg = float(0.5 * np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), K), P), y))

        # Ige = 1/2*y'*P*K*P*P*y;
        Ige = float(0.5 * np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), P), y))

        # Iee = 1/2*y'*P*P*P*y;
        Iee = (0.5 * np.dot(np.dot(np.dot(np.dot(y.T, P), P), P), y))

    elif alg == 1:
        # expected information
        # Igg = 1/2*trace(P*K*P*K)
        Igg = float(0.5 * np.trace(np.dot(np.dot(np.dot(P, K), P), K)))

        # Ige = 1/2*trace(P*K*P)
        Ige = float(0.5 * np.trace(np.dot(np.dot(P, K), P)))

        # Iee = 1/2*trace(P*P)
        Iee = (0.5 * np.trace(np.dot(P, P)))

    elif alg == 2:
        # observed information
        # Igg = -1/2*trace(P*K*P*K)+y'*P*K*P*K*P*y
        Igg = -float(0.5 * np.trace(np.dot(np.dot(np.dot(P, K), P), K))) + np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), K), P), y)

        # Ige = -1/2*trace(P*K*P)+y'*P*K*P*P*y
        Ige = -float(0.5 * np.trace(np.dot(np.dot(P, K), P))) + np.dot(np.dot(np.dot(np.dot(np.dot(y.T, P), K), P), P), y)

        # Iee = -1/2*trace(P*P)+y'*P*P*P*y;
        Iee = -(0.5 * np.trace(np.dot(P, P))) + np.dot(np.dot(np.dot(np.dot(y.T, P), P), P), y)

    # construct the score vector and the information matrix
    I = np.matrix([[Igg, Ige], [Ige, Iee]])

    # inverse of the information matrix
    invI = la.inv(I)
    # try:
    #     invI = la.inv(I)
    # except np.linalg.LinAlgError as e:
    #     invI = la.pinv(I)

    # standard error estimate
    # SE = sqrt((m2/Va)^2*((1-m2)^2*invI(1,1)-2*(1-m2)*m2*invI(1,2)+m2^2*invI(2,2)));
    SE = np.sqrt(
        ((m2 / Va)**2) *
        (
            (((1 - m2)**2) * invI[0, 0]) -
            (2 * (1 - m2) * m2) * invI[0, 1] +
            (m2**2) * invI[1, 1])
    )

    # --- diagnosis of convergence
    if ((iter == MaxIter) and (abs(Lnew - Lold) >= tol)):
        flag = 0
    else:
        flag = 1

    return {'flag': flag, 'm2': m2, 'SE': SE, 'Va': Va, 'Ve': Ve, 'Lnew': Lnew, 'dist': abs(Lnew - Lold)}


def variance_component_model(y, X, K, alg=0, tol=1e-4, MaxIter=100, verbal=False):
    '''
    Compute (using variance components moddel) the variation of y explained by
    the covariance matrix K, after regressing the effects of confounds in X.

    Args:
        y (ndarray):
            Nsubj x 1 vector of phenotypes (trait values)
        X (ndarray):
            Nsubj x Ncov covariate matrix (that contains "nuisance variables
            such as age, sex, site dummy, etc)
        K (ndarray):
            Nsubj x Nsubj anatomical similarity matrix (ASM)
            Note: K has to be a symmetric, positive semi-definite matrix with
            its diagonal elements averaging to 1. If K is not non-negative
            definite, its negative eigenvalues will be set to zero, and a warning
            will be printed.
        alg (int):
            the algorithm for the ReML; default=0
                alg = 0 - use the average information
                alg = 1 - use the expected information (Fisher scoring)
                alg = 2 - use the observed information
        tol (float):
            the tolerance for the convergence of the ReML algorithm; default=1e-4
        MaxIter (int):
            the maximum number of iterations for the ReML algorithm; default=100
        verbal (bool):
            a flag to set the verbosity; default=False

    Returns:
        dict:
            {
                'flag': flag indicates the convergence of the ReML algorithm (1 if
                        converged, and 0 if not),
                'm2': the morphometricity estimate (variance explained value),
                'SE': the standard error of the morphometricity estimate,
                'Va': the total anatomical/morphological variance,
                'Ve': the residual variance,
                'Lnew': the ReML likelihood when the algorithm is converged,
                'dist': difference of the last two steps in ReML convergence,
            }
    '''
    try:
        return(_variance_component_model_reimplemented(y, X, K, alg, tol, MaxIter, verbal))
    except Exception as e:
        print(e)
        return {'flag': 0, 'm2': np.nan, 'SE': np.nan, 'Va': np.nan, 'Ve': np.nan, 'Lnew': np.nan, 'dist': np.nan}
