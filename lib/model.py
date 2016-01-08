import random
import math, sys
import numpy as np
from scipy import optimize
from scipy.optimize import fminbound
import scipy.stats
from scipy.stats import kstest
from scipy.stats import chisquare

class T2_SSPSC:
    """
    This is the random variable that models the coalescence time of two different genes for the case
    of a instantaneous bottleneck. In this case, the bottleneck occurred at time T and the population size was
    multiplied by alpha. If not specified, the defaults values are aplha=0.1, T=1
    """
    def __init__(self, alpha=0.1, T=1):
        self.alpha= alpha
        self.T= T
    
    def pdf(self, t):
        """
        This is the density function of the random variable. The function receives a value t
        and returns the density function on this value
        """
        [alpha, T] = [self.alpha, self.T]
        t = np.array(t)
        return (t>0)*(t<T)*np.exp(-t) + (t>0)*(t>=T)*(1.0/alpha)*np.exp(-T- np.float_(t-T)/alpha)
    
    def cdf(self,t):
        """
        This function computes the Cumulative Distribution Function of the random variable
        """
        [alpha, T] = [self.alpha, self.T]
        t = np.array(t)
        return (t>0)*(t<T)*(1 - np.exp(-t)) + (t>0)*(t>=T)*(1 - np.exp(-T- np.float_(t-T)/alpha))
        
    def simulate_values(self, nb_of_observations):
        """
        This method simulates n independent values of the random variable
        """
        obs = []
        alpha=self.alpha
        T=self.T
        for i in range(nb_of_observations):
            X1 = random.expovariate(1)
            if X1<T:
                obs.append(X1)
            else:
                X2 = T + random.expovariate(1.0/alpha)
                obs.append(X2)
        return obs
  
    def __gi(self, ti, alpha, T):
        """
        This is the density function of the random variable evaluated at "ti". 
        """
        if (0<=ti<T):
            return math.exp(-ti)
        else:
            return math.exp(-T-float(ti-T)/alpha)/alpha

    def __log_gi(self, ti, alpha, T):
        """
        This is the log of the density function of the random variable evaluated at "ti".
        """
        if (0<=ti<T):
            return -ti
        else:
            return math.log(1.0/alpha) -T - float(ti-T)/alpha
 
    def log_likelihood(self, obs, alpha='', T=''):
        """
        This function returns the likelihood of obs. Here obs is a vector of observations
        """
        if alpha == '':
            alpha = self.alpha
        if T == '':
            T = self.T
        if (alpha<=0) or (T<0):
            return -sys.maxint
        temp = 0
        for ti in obs:
            temp+=self.__log_gi(ti, alpha, T)
        return temp

    def fast_log_likelihood(self, sorted_obs, alpha, T):
        """
        This function returns the likelihood of the observations. It assumes that the vector of observations
        is sorted descent
        """
        if (alpha<=0) or (T<0):
            return -sys.maxint
        temp=0
        for i in range(len(sorted_obs)):
            if sorted_obs[i]>=T:
                temp+=math.log(1.0/alpha) -T - float(sorted_obs[i]-T)/alpha
            else: break
        if sorted_obs[i]>=T: # This is for avoiding to substract the last i when all are >= T
            i+=1
        return temp-sum(sorted_obs[i:])

    def right_limit_fast_log_likelihood(self, sorted_obs, alpha, T):
        """
        The likelihood function has a jump discontinuity at every ti. It is continuous from the left but 
        discontinuous from the right. This function computes the right limit of the likelihood function.
        """
        if (alpha<=0) or (T<0):
            return -sys.maxint
        temp=0
        for i in range(len(sorted_obs)):
            if sorted_obs[i]>T:
                temp+=math.log(1.0/alpha) -T - float(sorted_obs[i]-T)/alpha
            else: break
        if sorted_obs[i]>T: # This is for avoiding to substract the last i when all are > T
            i+=1
        return temp-sum(sorted_obs[i:])
    
    def max_likelihood_estimation(self, obs, x0=''):
        """
        This method estimates the parameters that maximize the likelihood. It is necessary to enter
        some initial guess x0. If the initial guess is not given, it takes alpha=1 and T=m, where 
        "m" is the mean of the observed values.
        """
        initial_guess=x0
        if initial_guess=='':
            initial_guess=np.array([1, sum(obs)/len(obs)])
        llB1 = lambda x: -self.log_likelihood(obs, x[0], x[1])
        res = optimize.fmin(llB1, initial_guess, disp=False)
        return [res[0], res[1], self.log_likelihood(obs, res[0], res[1])]
    
    def compute_alpha(self, obs, t):
        """
        This method computes the value of alpha for the exact computation
        of the likelihood. 
        obs should be a decreasing ordered list and at least
        obs[0]>t
        """
        n_temp = 0
        sum_temp = 0
        obs.sort(reverse=True)
        for i in range(len(obs)):
            if obs[i] < t: break
            n_temp += 1
            sum_temp += obs[i]
        return float(sum_temp - n_temp*t)/(n_temp)

    def compute_right_alpha(self, obs, t):
        """
        This method computes the value of alpha for the exact computation
        of the likelihood when using the right limit formulae 
        obs should be a decreasing ordered list and at least
        obs[0]>t
        """
        n_temp = 0
        sum_temp = 0
        obs.sort(reverse=True)
        for i in range(len(obs)):
            if obs[i] <= t: break
            n_temp += 1
            sum_temp += obs[i]
        return float(sum_temp - n_temp*t)/(n_temp)

    def exact_maxllk(self, obs):
        # This function uses a exact method to find the maximum of the log_likelihood function
        obs.sort(reverse=True)
        # Computing the first value
        best_t = obs[0]
        best_alpha = 1 # The likelihood function doesn't depend on alpha beyond the biggest T
        sup_llk = self.fast_log_likelihood(obs, best_alpha, best_t)
        best_llk = sup_llk # The function is continuous in T = max(obs)
        for ti in obs[1:]+[0]:
            alpha_i_left = self.compute_alpha(obs, ti)
            alpha_i_right = self.compute_right_alpha(obs, ti)
            llk_i_left = self.fast_log_likelihood(obs, alpha_i_left, ti)
            llk_i_right = self.right_limit_fast_log_likelihood(obs, alpha_i_right, ti)
            if llk_i_left > best_llk:
                best_t = ti
                best_alpha = alpha_i_left
                best_llk = llk_i_left
            if llk_i_right > best_llk:
                best_t = ti
                best_alpha = alpha_i_right
                best_llk = llk_i_right
        return [best_alpha, best_t, best_llk]
        
    def exact_maxllk_integer(self, obs):
        """
        Do the same algorithm that 'exact_maxllk' but considering that alpha
        is integer. The algorithm is roughly the same due to the form of the 
        likelihood function (i.e. for constant values of T, the function 
        of alpha is convex)
        """
        obs.sort(reverse=True)
        # Computing the first value
        best_t = obs[0]
        best_alpha = 1 # The likelihood function doesn't depend on alpha beyond the biggest T
        sup_llk = self.fast_log_likelihood(obs, best_alpha, best_t)
        best_llk = sup_llk # The function is continuous in T = max(obs)
        for ti in obs[1:]+[0]:
            # The maximum of a function is over one of the observed values of 
            # coalescence times. 
            # Moreover, for constant values of T, the likelihood function in 
            # convex. We get then the better alpha that is integer
            alpha_i_left_real = self.compute_alpha(obs, ti)
            (alpha_inf, ignore) = divmod(alpha_i_left_real, 1)
            alpha_sup = alpha_inf + 1
            llk_inf = self.fast_log_likelihood(obs, alpha_inf, ti)
            llk_sup = self.fast_log_likelihood(obs, alpha_sup, ti)
            if llk_inf > llk_sup:
                alpha_i_left = alpha_inf
            else:
                alpha_i_left = alpha_sup
            
            alpha_i_right_real = self.compute_right_alpha(obs, ti)
            (alpha_inf, ignore) = divmod(alpha_i_right_real, 1)
            (alpha_inf, ignore) = divmod(alpha_i_right_real, 1)
            alpha_sup = alpha_inf + 1
            llk_inf = self.fast_log_likelihood(obs, alpha_inf, ti)
            llk_sup = self.fast_log_likelihood(obs, alpha_sup, ti)
            if llk_inf > llk_sup:
                alpha_i_right = alpha_inf
            else:
                alpha_i_right = alpha_sup
            
            llk_i_left = self.fast_log_likelihood(obs, alpha_i_left, ti)
            llk_i_right = self.right_limit_fast_log_likelihood(obs, alpha_i_right, ti)
            if llk_i_left > best_llk:
                best_t = ti
                best_alpha = alpha_i_left
                best_llk = llk_i_left
            if llk_i_right > best_llk:
                best_t = ti
                best_alpha = alpha_i_right
                best_llk = llk_i_right
        return [best_alpha, best_t, best_llk]        

class T2_StSI:
    """
    This is the random variable that models the coalescence time of two different genes for the case
    of a structured population with n islands and migration rate M. If the parameters are not set the 
    default values are n=10 and M=0.1
    The values A, B, and E are computed from n and M when an object is created.
    """  
    def __init__(self, n=10, M=0.1):
        self.__n = n
        self.__M = M
        self.__update_constants()    
    
    def pdf(self, t):
        """
        This is the density function of the random variable
        """
        [a, alpha, beta] = [self.__a, self.__alpha, self.__beta]
        return (t>0)*(a*np.exp(-alpha*t) + (1-a)*np.exp(-beta*t))
    
    def cdf(self, t):
        """
        This function computes the cumulative distribution function of the random variable
        """
        [a, alpha, beta] = [self.__a, self.__alpha, self.__beta]
        return (t>0)*(1 - a*np.exp(-alpha*t)/alpha - (1-a)*np.exp(-beta*t)/beta)

    def __update_constants(self):
        """
        This method computes the constants (intermediate values) used for computations in
        simulating values and for computing the likelihood.
        """
        [self.__A, self.__B, self.__a, self.__alpha, self.__beta, self.__E, self.__AplusB, self.__AminusB]=self.__compute_constants(self.__n,self.__M)

    def __compute_constants(self, n, M):
        A = 1+(float(n)/(n-1))*M
        B = math.sqrt(A**2-float(4*M)/(n-1))
        a = 0.5 + (1.0 + (float(n-2)*M)/(n-1))/(2*B)
        alpha = 0.5*(A+B)
        beta = 0.5*(A-B)
        E = (1 + (float(n-2)*M)/(n-1)) / (2*B)
        AplusB = A+B
        AminusB = A-B
        return [A, B, a, alpha, beta, E, AplusB, AminusB]
        
    def compute_constants(self, n, M):
        A = 1+(float(n)/(n-1))*M
        B = math.sqrt(A**2-float(4*M)/(n-1))
        a = 0.5 + (1.0 + (float(n-2)*M)/(n-1))/(2*B)
        alpha = 0.5*(A+B)
        beta = 0.5*(A-B)
        E = (1 + (float(n-2)*M)/(n-1)) / (2*B)
        AplusB = A+B
        AminusB = A-B
        return [A, B, a, alpha, beta, E, AplusB, AminusB]

    def __compute_gi(self, obs_i, n, M):
        """
        Compute some intermediate values for every observation.
        This values are used for the likelihood computation
        """
        if (n<>self.__n) or M<>self.__M:
            [A, B, a, alpha, beta, E, AplusB, AminusB] = self.__compute_constants(n, M)
        else:
            AplusB=self.__AplusB
            AminusB=self.__AminusB
            E = self.__E
        Ci = math.exp(-float(obs_i)*AplusB/2) + math.exp(-float(obs_i)*AminusB/2)
        Di = math.exp(-float(obs_i)*AplusB/2) - math.exp(-float(obs_i)*AminusB/2)
        return 0.5*Ci + E*Di
        
    # We use properties for changing the values of n and M, given that changing those
    # values implies changing other constants used in computations
    def get_n(self):
        return self.__n
    def set_n(self, n):
        self.__n = n
        self.__update_constants()
    n = property(get_n, set_n)
    
    def get_M(self):
        return self.__M
    def set_M(self, M):
        self.__M = M
        self.__update_constants()
    M = property(get_M, set_M)
        
    def simulate_values(self, nb_of_observations):
        """
        Simulate n independent values of the random variable
        """
        observations = []
        for i in range(nb_of_observations):
            if random.random() < (float(self.__a)/self.__alpha):
                observations.append(random.expovariate(self.__alpha))
            else:
                observations.append(random.expovariate(self.__beta))
        return observations  
        
    def log_likelihood(self, obs, n='', M=''):
        """
        Compute the likelihood of the observations given the parameters
        n and M. If n and M are not entered it uses the parameters of the actual instance. 
        """
        if n=='':
            n=self.n
        if M=='':
            M=self.M
        temp = 0
        if (n<2) or (M<=0):
            return -sys.maxint
        #terms = []
        for ti in obs:
            # Here, for ti very big, gi becomes 0. That's a problem for the logarithm
            gi = self.__compute_gi(ti, n, M)
            try:
                temp += math.log(gi)
            except:
                temp -=ti
            #terms.append(math.log(a*exp(-alpha*t)+(1-a)*exp(-beta*t)))
        return temp
    
    def __bound_log_likelihood(self, obs, n, M, bound, direction):
        """
        This function allows to compute the likelihood of the parametes in the case where it is necessary to
        define some bounds for the value of n
        """
        if (direction == 'left') and (n>bound):
            return -sys.maxint
        elif (direction == 'right') and (n<bound):
            return -sys.maxint
        return self.log_likelihood(obs, n, M)
    
    def bound_log_likelihood(self, obs, n, M, bound, direction):
        """
        This function allows to compute the likelihood of the parametes in the case where it is necessary to
        define some bounds for the value of n
        """
        if (direction == 'left') and (n>bound):
            return -sys.maxint
        elif (direction == 'right') and (n<bound):
            return -sys.maxint
        return self.log_likelihood(obs, n, M)
    
    def max_likelihood_estimation(self, obs, min_M = 0.001, max_M = 1000):
        # This method finds an aproximation to the Maximum of the
        # likelihood function. The algorithm is too simple and
        # may give a very bad solution for some cases.
        # The strategy is: Just find a maximum for the function (this
        # give two real numbers. Then we round the "n" parameters and we
        # return the bigger evaluation of f(n, M) and f(n+1, M)
        ll = lambda x: -self.log_likelihood(obs, x[0], x[1])
        initial_guess = np.array([2, 1])
        res_temp = optimize.fmin(ll, initial_guess, disp=False)
        solution_temp = res_temp
        if np.trunc(solution_temp[0]) == solution_temp[0]:
            return [solution_temp[0], solution_temp[1], self.log_likelihood(obs, solution_temp[0], solution_temp[1])]
        else:
            M_found = solution_temp[1]
            n_left = np.trunc(solution_temp[0])
            n_right = n_left+1
            sol_left = self.log_likelihood(obs, n_left, M_found)
            sol_right = self.log_likelihood(obs, n_right, M_found)
            if sol_left > sol_right:
                return [n_left, M_found, sol_left]
            else:
                return [n_right, M_found, sol_right]
                
class T2_StSId(T2_StSI):
    def __init__(self, n=10, M=0.1):
        T2_StSI.__init__(self, n, M)
        [A, B, a, alpha, beta, E, AplusB, AminusB] = \
        T2_StSI.compute_constants(self, n, M)
        [self.a, self.alpha, self.beta] = [a, alpha, beta]
        self.gamma = np.true_divide(M, n-1)
        self.c = np.true_divide(self.gamma, beta-alpha)
    
    def pdf(self, t):
        [alpha, beta, c] = [self.alpha, self.beta, self.c]
        return (t>0)*(np.exp(-alpha*t) - np.exp(-beta*t))*c
    
    def cdf(self, t):
        [alpha, beta] = [self.alpha, self.beta]
        return (t>0)*(1 - np.true_divide(alpha*np.exp(-beta*t) - 
                        beta*np.exp(-alpha*t), alpha-beta))
    
    def simulate_values(self, nb_of_observations):
        print("not implemented")
        return False

# The two next clases are for comparing the simulated values produced by
# T2_SSPSC and T2_StSI with the equivalent ms commands. Here we have to take
# into accound that ms times are scaled by 4N and our times are scaled by
# 2N (N is the population size at the present)

class T2_SSPSC_MS(T2_SSPSC):
    """
    This class should produce the same values that the command
    ms 2 200000 -T -L -eN T alpha
    where "alpha" and "T" are the parameters of the T2_SSPSC
    distribution function 
    """
    def __init__(self, alpha, T):
        T2_SSPSC.__init__(self, alpha, T*2)
    def cdf(self, t):
        return super(T2_SSPSC_MS, self).cdf(t*2)
    def pdf(self, t):
        return super(T2_SSPSC_MS, self).pdf(t*2)
    def simulate_values(self, nb_of_observations):
        """
        This method simulates n independent values of the random variable
        We are going to divide by two the rates of the exponential distributions
        because we are using a scale two times bigger (T = 4N) as in MS.
        """
        temp = super(T2_SSPSC_MS, self).simulate_values(nb_of_observations)
        return np.true_divide(np.array(temp),2)

class T2_StSI_MS(T2_StSI):
    """
    This class should produce the same values that the command
    ms 2 20000 -T -L -I 9 2 0 0 0 0 0 0 0 0 0.1
    for n=9 islands and M=4Nm=0.1
    As in our function M=2Nm, we need to initialize our distribution
    with M=0.1/2=0.05
    """
    def __init__(self, n, M):
        T2_StSI.__init__(self, n, M)
    def cdf(self, t):
        return super(T2_StSI_MS, self).cdf(t*2)
    def pdf(self, t):
        return super(T2_StSI_MS, self).pdf(t*2)
    def simulate_values(self, nb_of_observations):
        temp = super(T2_StSI_MS, self).simulate_values(nb_of_observations)
        return np.true_divide(np.array(temp), 2)

class Comparator:
    """
    This class is for comparing the results of the simulations.
    """
    def Akaike_IC(self, lnl, k):
        """
        Akaike information criterion. lnl is the log of the likelihood and k is the number of parameters
        of the model
        """
        return 2 * k - 2 * lnl
    
    def AIC_compare(self, max_log_Lik_T2_SSPSC, max_log_Lik_T2_StSI):
        """
        Compares both models following AIC. 
        Returns a tuple in the form (0,v) or (1,v).
        0 - the best model is T2_SSPSC
        1 - the best model is T2_StSI
        the second value of the tuple is the relative probability that the other method minimizes the estimated
        information loss.
        Ex: (0, 0.2) means that the selected model was T2_SSPSC and that T2_StSI has probability of 0.2
        of estimate the information loss
        """
        T2_SSPSC_AIC = self.Akaike_IC(max_log_Lik_T2_SSPSC,2)
        T2_StSI_AIC = self.Akaike_IC(max_log_Lik_T2_StSI, 2)
        if T2_SSPSC_AIC < T2_StSI_AIC:
            best_model = 0
            relative_prob = np.exp(np.true_divide(T2_SSPSC_AIC-T2_StSI_AIC,2))
        else:
            best_model = 1
            relative_prob = np.exp(np.true_divide(T2_StSI_AIC-T2_SSPSC_AIC,2))
        return (best_model, relative_prob)

    def single_exp(self, options):
        """
        This method returns a line in the form 
        Original_variable | real_parameters | number_of_observations | log-likelihood_of_real_params | Estim_params_T2_SSPSC | log-likelihood_T2_SSPSC | p-value_T2_SSPSC | Estim_params_T2_StSI | log-likelihood_T2_StSI | p-value_T2_StSI | AIC_selected_model | AIC_relative_prob
        The input is an array of options containing:
        ['choice of distrib', 'param1, param2', 'number of values']
        The input is going to be used to know the kind of experiment
        we are going to do. 
        The procedure is to simulate data with one distribution and then
        to do a KS test for determining if the data comes from a
        SSPSC or from a StSI model.
        Moreover, it computes AIC values in order to say wich model explains better the simulated data
        """
        params = options[1].split(',')
        if options[0]=='1':
            alpha=float(params[0])
            T=float(params[1])
            X = T2_SSPSC(alpha, T)
            type_of_variable = 'T2_SSPSC'
            real_params = '({}, {})'.format(alpha, T)
        elif options[0]=='2':
            n=int(params[0])
            M=float(params[1])
            X = T2_StSI(n, M)
            type_of_variable = 'T2_StSI'
            real_params = (n, M).__str__()
        number_of_observations = int(options[2])
        alpha_integer = options[3]
        
        # We are going to simulate "number_of_observations" independent values. Then we use a half of that values for parameters estimation and the other half for the KS-test
        
        # Simulating values
        obs_estim = X.simulate_values(number_of_observations/2)
        obs_test = X.simulate_values(number_of_observations/2)
        obs = obs_estim + obs_test
        real_likelihood = X.log_likelihood(obs_estim)

        # Estimating parameters of both models and doing a KS test
        B = T2_SSPSC() # We initialize a variable with default parameters
        if int(alpha_integer) == 0:
            [alpha, T, ll_fittedB] = B.exact_maxllk(obs_estim)
        else:
            [alpha, T, ll_fittedB] = B.exact_maxllk_integer(obs_estim)
        fittedB = T2_SSPSC(alpha, T)
        KS_fittedB = kstest(obs_test, fittedB.cdf)
        pvalue_fittedB = KS_fittedB[1]
        S = T2_StSI() # We initialize a variable with default parameters
        [n, M, ll_fittedS] = S.max_likelihood_estimation(obs_estim)
        fittedS = T2_StSI(n, M)
        KS_fittedS = kstest(obs_test, fittedS.cdf)
        pvalue_fittedS = KS_fittedS[1]

        # Computing AIC values
        (best_model, relative_prob) = self.AIC_compare(ll_fittedB, ll_fittedS)

        result_text= type_of_variable+' | '+real_params+' | '+options[2]+' | '+real_likelihood.__str__()+' | '+(fittedB.alpha, fittedB.T).__str__()+' | '+ll_fittedB.__str__()+' |  '+pvalue_fittedB.__str__()+' | '+(fittedS.n, fittedS.M).__str__()+' | '+ll_fittedS.__str__()+' | '+pvalue_fittedS.__str__()+' | '+best_model.__str__()+' | '+relative_prob.__str__()
        observations = obs.__str__()
        return [result_text, observations]

    def multiple_experiments(self, options):
        '''
        This method returns the result of doing multiple single_experiments. It calls the method 'single_experiment'
        with the necessary parameters.
        The input of this method is an array of string in the form:
        ['choice of distrib', 'list of params1', 'list of params2', 'list of number of values']
        '''
        output_text = 'Original_variable | real_parameters | number_of_observations | log-likelihood_of_real_params | Estim_params_T2_SSPSC | log-likelihood_T2_SSPSC | p-value_T2_SSPSC | Estim_params_T2_StSI | log-likelihood_T2_StSI | p-value_T2_StSI | AIC_selected_model | AIC_relative_prob\n'
        observations_text = ''
        single_exp_options = [options[0], '', '', options[4]]
        list_param1 = options[1].split(',')
        list_param2 = options[2].split(',')
        list_number_of_values = options[3].split(',')
        for number_of_values in list_number_of_values:
            single_exp_options[2]=number_of_values
            temp = ''
            for param1 in list_param1:
                for param2 in list_param2:
                    single_exp_options[1] = '{},{}'.format(param1, param2)
                    [result, observations] = self.single_exp(single_exp_options)
                    output_text += result+'\n'
                    observations_text += observations+'\n'
        return [output_text, observations_text]

class Parser:
    """
    This class is used to parse the configuration file and return the parameters
    """
    def __init__(self, filename):
        self.filename = filename
    
    def get_array_of_parameters(self):
        a = open(self.filename)
        text = a.read()
        a.close()
        params = text.split('\n')
        options = []
        for i in params:
            if (len(i)>0) and (i[0]<>'#'):
                options.append(i)
        # Reading options
        return options

