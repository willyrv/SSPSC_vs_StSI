

class Experiment:
    def __init__(self, Orig_dist, param_1, param_2, numb_of_obs, llk_rp, T2_SSPSC_est_alpha, T2_SSPSC_est_T, llk_T2_SSPSC, T2_SSPSC_pvalue, T2_StSI_est_n, T2_StSI_est_M, llk_T2_StSI, T2_StSI_pvalue, AIC_choice, AIC_rel_prob, rejection_bound=0.05):
        self.Orig_dist = Orig_dist
        self.param_1 = param_1
        self.param_2 = param_2
        self.numb_of_obs = numb_of_obs
        self.llk_rp = llk_rp
        self.T2_SSPSC_est_alpha = T2_SSPSC_est_alpha
        self.T2_SSPSC_est_T = T2_SSPSC_est_T
        self.llk_T2_SSPSC = llk_T2_SSPSC
        #self.sup_llk_T2_SSPSC = sup_llk_T2_SSPSC
        self.T2_SSPSC_pvalue = T2_SSPSC_pvalue
        self.T2_StSI_est_n = T2_StSI_est_n
        self.T2_StSI_est_M = T2_StSI_est_M
        self.llk_T2_StSI = llk_T2_StSI
        self.T2_StSI_pvalue = T2_StSI_pvalue
        self.orig_dist_code = ''
        if (Orig_dist == "T2_StSI") or (Orig_dist.__contains__("Structure")):
            self.orig_dist_code=1
            self.n = param_1
            self.M = param_2
        if (Orig_dist == "T2_SSPSC") or (Orig_dist.__contains__("Bottleneck")):
            self.orig_dist_code=0
            self.alpha = param_1
            self.T = param_2
        self.AIC_choice = AIC_choice
        self.AIC_rel_prob = AIC_rel_prob
        
        # Results of the KS test
        # list all possible cases
        self.KS_0rej = False
        self.KS_2rej = False
        self.KS_select = False
        #self.KS_right_choice = False
        #self.KS_wrong_choice = False
        
        if self.T2_SSPSC_pvalue<rejection_bound:
            if self.T2_StSI_pvalue<rejection_bound:
                # KS rejects both models
                self.KS_2rej = True
            else:
                # KS selects T2_StSI
                self.KS_select = True
                self.KS_right_choice = (self.orig_dist_code == 1)
                self.KS_wrong_choice = not self.KS_right_choice
        else:
            if self.T2_StSI_pvalue<rejection_bound:
                # KS selects T2_SSPSC
                self.KS_select = True
                self.KS_right_choice = (self.orig_dist_code == 0)
                self.KS_wrong_choice = not self.KS_right_choice
            else:
                # KS doesn't reject any model
                self.KS_0rej = True
                
        # Results of the AIC
        self.AIC_right_choice = (self.orig_dist_code == self.AIC_choice)
        self.AIC_wrong_choice = not self.AIC_right_choice

    def reject_T2_SSPSC(self, alpha):
        if self.T2_SSPSC_pvalue < alpha:
            return True
        else:
            return False

    def reject_T2_StSI(self, alpha):
        if self.T2_StSI_pvalue < alpha:
            return True
        else:
            return False

class FullExperiment:
    """
    This class will contain the results of the experiments in a list of objetcs of type "Experiment"
    """
    def __init__(self,configfile, experiment_results_filename):
        self.configfile = configfile
        # We read all the parameters from the config file and stock them in class variables
        self.Configure(configfile)
        self.header_line = "EXP_No | Type_of_variable | real_parameters | number_of_observations | log-likelihood_of_real_params | Params_T2_SSPSC | log-likelihood_T2_SSPSC | sup_T2_SSPSC_likelihood | p-value_T2_SSPSC | params_T2_StSI | log-likelihood_T2_StSI | p-value_T2_StSI"
        self.load_Experiments(experiment_results_filename)
    
    def Configure(self, configfile):
        # This method can be used to change the config file and all the configuration parameters
        self.configfile = configfile
        a = open(configfile)
        text = a.read()
        a.close()
        params = text.split('\n')
        options = []
        for i in params:
            if (len(i)>0) and (i[0]<>'#'):
                options.append(i)
        # Reading options
        if options[0]==1:
            self.original_dist='T2_SSPSC'
        elif options[0]==2:
            self.original_dist='T2_StSI'
        self.array_of_number_obs = options[3]
        self.array_of_param1 = options[1]
        self.array_of_param2 = options[2]
        self.list_of_exp = []        

    def add_Experiment(self, line_result, obs=''):
        l = line_result.split('|')
        Orig_dist = l[0].strip()
        p = l[1].strip()
        p = p.strip('()')
        p = p.split(',')
        param_1 = float(p[0])
        param_2 = float(p[1])
        numb_of_obs = int(l[2].strip())
        llk_rp = float(l[3])
        p = l[4].strip()
        p = p.strip('()')
        p = p.split(',')
        T2_SSPSC_est_alpha = float(p[0])
        T2_SSPSC_est_T = float(p[1])
        llk_T2_SSPSC = float(l[5].strip())
        #sup_llk_T2_SSPSC = float(l[7].strip())
        T2_SSPSC_pvalue = float(l[6].strip())
        p = l[7].strip()
        p = p.strip('()')
        p = p.split(',')
        T2_StSI_est_n = float(p[0])
        T2_StSI_est_M = float(p[1])
        llk_T2_StSI = float(l[8].strip())
        T2_StSI_pvalue = float(l[9].strip())
        AIC_choice = int(l[10])
        AIC_rel_prob = float(l[11])
        self.list_of_exp.append(Experiment(Orig_dist, param_1, param_2, numb_of_obs, llk_rp, T2_SSPSC_est_alpha, T2_SSPSC_est_T, llk_T2_SSPSC, T2_SSPSC_pvalue, T2_StSI_est_n, T2_StSI_est_M, llk_T2_StSI, T2_StSI_pvalue, AIC_choice, AIC_rel_prob))

    def load_Experiments(self, filename):
        # This method loads the results of the experiments from a file
        a = open(filename, 'r')
        l = a.readline() # This is for reading the first line
        l = a.readline()
        while l<>'\n':
            self.add_Experiment(l)
            l = a.readline()
        a.close()
        
    def save2file(self, filename):
        a = open(filename, 'w')
        pickle.dump(self, a)
        a.close()
        
    def KS_resume(self, short_resume=False):
        """
        Prints the number (and the percent) of times the KS has done a good choice as well
        as the proportion of no choices
        """
        temporal_text = ['Total | nL | KS_choices (total, percent) | KS_right | KS_wrong | KS_0rej | KS_2rej']
        n_obs = self.array_of_number_obs.split(',')
        n_obs = [int(i) for i in n_obs]
        for nL in n_obs:
            Total = 0
            KS_choices = 0
            KS_right = 0
            KS_wrong = 0
            KS_0rej = 0
            KS_2rej = 0
            for exp in self.list_of_exp:
                if exp.numb_of_obs==nL:
                    Total+=1
                    if exp.KS_select:
                        KS_choices+=1
                        KS_right+=exp.KS_right_choice
                        KS_wrong+=exp.KS_wrong_choice
                    KS_0rej+=exp.KS_0rej
                    KS_2rej+=exp.KS_2rej
            percent_KS_choices = float(KS_choices)/Total * 100
            percent_KS_right = float(KS_right)/Total * 100
            rel_percent_KS_right = float(KS_right)/KS_choices * 100
            percent_KS_wrong = float(KS_wrong)/Total * 100
            rel_percent_KS_wrong = float(KS_wrong)/KS_choices * 100
            percent_KS_0rej = float(KS_0rej)/Total * 100
            percent_KS_2rej = float(KS_2rej)/Total * 100
            
            if short_resume:
                # do short 
                not_implemented = 1
            else:
                temporal_text.append('{} | {} | {}:{}% | {}:{}:{}% | {}:{}:{}% | {}:{}% | {}:{}%'.format(
                                    Total, nL, KS_choices, percent_KS_choices, KS_right, percent_KS_right, 
                                    rel_percent_KS_right, KS_wrong, percent_KS_wrong, rel_percent_KS_wrong, 
                                    KS_0rej, percent_KS_0rej, KS_2rej, percent_KS_2rej))
        output_text = '\n'.join(temporal_text)
        print output_text
    
    def AIC_resume(self):
        """
        Computes the number of times the AIC does the right choice when KS fails
        """
        temporal_text = ['Total | nL | KS_0rej (total:percent) | KS_2rej | AIC_right_choice (total:relative_percent:absolute_percent) | AIC_wrong_choice']
        n_obs = self.array_of_number_obs.split(',')
        n_obs = [int(i) for i in n_obs]
        for nL in n_obs:
            Total = 0
            KS_0rej = 0
            KS_2rej = 0
            AIC_right_choice = 0
            AIC_wrong_choice = 0
            for exp in self.list_of_exp:
                if exp.numb_of_obs==nL:
                    Total+=1
                    if not exp.KS_select:
                        KS_0rej+=exp.KS_0rej
                        KS_2rej+=exp.KS_2rej
                        AIC_right_choice+=exp.AIC_right_choice
                        AIC_wrong_choice+=exp.AIC_wrong_choice

            percent_KS_0rej = float(KS_0rej)/Total * 100
            percent_KS_2rej = float(KS_2rej)/Total * 100
            rel_percent_AIC_right_choice = float(AIC_right_choice)/(KS_0rej+KS_2rej) * 100
            percent_AIC_right_choice = float(AIC_right_choice)/Total * 100
            rel_percent_AIC_wrong_choice = float(AIC_wrong_choice)/(KS_0rej+KS_2rej) * 100
            percent_AIC_wrong_choice = float(AIC_wrong_choice)/Total * 100
            

            temporal_text.append('{} | {} | {}:{}% | {}:{}% | {}:{}:{}% | {}:{}:{}%'.format(
                Total, nL, KS_0rej, percent_KS_0rej, KS_2rej, percent_KS_2rej, 
                AIC_right_choice, rel_percent_AIC_right_choice, percent_AIC_right_choice, AIC_wrong_choice, rel_percent_AIC_wrong_choice, percent_AIC_wrong_choice))
        output_text = '\n'.join(temporal_text)
        print output_text

