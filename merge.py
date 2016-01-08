#!/usr/bin/env python

OUTPUT_FILENAME = "output"
TOTAL_OF_FILES = 2

if __name__ == "__main__":
    # make an array of files of size 100
    results = []
    for i in range(1, TOTAL_OF_FILES+1):
        a = open("./{}{}.txt".format(OUTPUT_FILENAME, i), 'r')
        text = a.read()
        a.close()
        text = text.split('\n')
        results.append(text[1:-2])
    text = "Original_variable | real_parameters | number_of_observations | log-likelihood_of_real_params | Estim_params_T2_SSPSC | log-likelihood_T2_SSPSC | p-value_T2_SSPSC | Estim_params_T2_StSI | log-likelihood_T2_StSI | p-value_T2_StSI | AIC_selected_model | AIC_relative_prob\n"
    for line in range(len(results[0])):
        for i in range(TOTAL_OF_FILES):
            text+=(i+1).__str__() + ' | ' + results[i][line]+'\n'
    a = open('100exp_OUTPUT', 'w')
    a.write(text)
    a.close()

    
