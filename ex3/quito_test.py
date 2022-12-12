# import sys
# import os
# import numpy
# curPath = os.path.abspath(os.path.dirname(__file__))
# rootPath = os.path.split(curPath)[0]
# sys.path.append(os.path.split(rootPath)[0])
#
# if sys.platform.startswith('win32'):
#     ROOT='\\'
# elif sys.platform.startswith('linux'):
#     ROOT='/'
# elif sys.platform.startswith('darwin'):
#     ROOT='/'
#
# K=200
# M=20
# BUDGET=0
# C_LEVEL=0.01
# S_TEST = 'one-sample Wilcoxon signed rank test'
# T1=0
# T2=0
# START = 0
#
# from qiskit import (
#     #IBMQ,
#     QuantumCircuit,
#     QuantumRegister,
#     ClassicalRegister,
#     execute,
#     Aer,
# )
# from rpy2 import robjects as robjects
# import numpy as np
# import os
#
# import configparser
# import os.path
# import importlib
# import time
# import warnings
#
# import logging
# from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
# rpy2_logger.setLevel(logging.ERROR)   # will display errors, but not warningswo
#
# def check_unique(l):
#     return len(l) == len(set(l))
#
# def end_running():
#     exit()
#
# def dec2bin(n,bit):
#     a = 1
#     list = []
#     while a > 0:
#         a, b = divmod(n, 2)
#         list.append(str(b))
#         n = a
#     s = ""
#     for i in range(len(list) - 1, -1, -1):
#         s += str(list[i])
#     s = s.zfill(bit)
#     return s
#
# def input_group(valid_input):
#     index = [] #unique input index
#     index_flag = valid_input[0]
#     index.append(0)
#     for i in range(1,len(valid_input)):
#         if valid_input[i] != index_flag:
#             index.append(i)
#             index_flag = valid_input[i]
#     return index
#
# def check_full_ps(valid_input, p):
#     index = input_group(valid_input)
#     p_sum = 0
#     for i in range(len(index)):
#         start = index[i]
#         if i == len(index) - 1:
#             end = len(valid_input)
#         else:
#             end = index[i+1]
#         for j in range(start,end):
#             p_sum += p[j]
#         if p_sum != 1:
#             print("Error: This is not a complete program specification.")
#             end_running()
#         else:
#             p_sum = 0
#
# def check_partial_ps(valid_input, valid_output, p):#check whether there is any input having full ps(sum of p of all outputs of one input is 1)
#     index = input_group(valid_input)
#     p_sum = 0
#     input_full = []
#     output_full = []
#     for i in range(len(index)):
#         start = index[i]
#         if i == len(index) - 1:
#             end = len(valid_input)
#         else:
#             end = index[i+1]
#         for j in range(start,end):
#             p_sum += p[j]
#         if p_sum == 1:
#             for k in range(start,end):
#                 input_full.append(valid_input[k])
#                 output_full.append(valid_output[k])
#         p_sum = 0
#     return input_full, output_full
#
# def get_unique(l):
#     unique = []
#     for i in l:
#         if i not in unique:
#             unique.append(i)
#     return unique
#
# def get_all(bit):
#     all = []
#     for i in range(pow(2,bit)):
#         i_bin = dec2bin(i, bit)
#         all.append(i_bin)
#     return all
#
# def execute_quantum_program(inputID, outputID, num_qubit, i, module_name):
#     q = QuantumRegister(num_qubit)
#     c = ClassicalRegister(len(outputID))
#     qc = QuantumCircuit(q, c)
#     for j in range(len(inputID)):
#         # print(i)
#         if i[len(inputID) - 1 - j] == '1':
#             # print(int(inputID[j]))
#             qc.x(int(inputID[j]))
#     module = importlib.import_module(module_name)
#     run_method = getattr(module,"run")
#     run_method(qc)
#     result = execute(qc, Aer.get_backend('qasm_simulator'), shots=1).result().get_counts(qc)
#     return result
#
# def check_same(l,value):
#     for i in range(len(l)):
#         if l[i] != value:
#             return False
#     return True
#
# def check_WOO(i, o , valid_inputs, valid_outputs):
#     flag = False
#     for k in range(len(valid_inputs)):
#         if valid_inputs[k] == i and valid_outputs[k] == o:
#             flag = True
#     if flag == False:
#         print('fail for woo')
#         return False
#     return True
#
# def wilcoxon(fre,p):
#     pvalue = []
#     for i in range(len(p)):
#         if np.isnan(fre[i]).any() == True:
#             pvalue.append(-1)
#         elif check_same(fre[i],p[i]) == True:
#             # print('exactly the same')
#             pvalue.append(1)
#         else:
#             fre_r = robjects.FloatVector(fre[i])
#             # print('fre='+str(fre[i]))
#             # print('p='+str(p[i]))
#             robjects.r('''
#                         wtest<-function(data,exp,c_level){
#                             test_result <- wilcox.test(data, mu = exp, conf.int = TRUE, conf.level = 1-c_level)
#                             pvalue = test_result$p.value
#                             return (pvalue)
#                         }
#             ''')
#             t = robjects.r['wtest'](fre_r, p[i], C_LEVEL)
#             pvalue.append(t[0])
#     #print(pvalue)
#     return pvalue
#
# def judge_ass_result(inputs, outputs, pvalue, f):
#     for i in range(len(pvalue)):
#         if pvalue[i] == -1:
#             f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + "Not enough inputs for statistical test.")
#             f.write('\n')
#         elif pvalue[i] < C_LEVEL:
#             f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + str(pvalue[i]) + "; ")
#             f.write("Result: Fail for OPO")
#             f.write('\n')
#         elif pvalue[i] >= C_LEVEL:
#             f.write("(" + str(inputs[i]) + "," + str(outputs[i]) + "): " + str(pvalue[i]) + "; ")
#             f.write("Result: Inconclusive")
#             f.write('\n')
#
# def process_bar(percent, start_str='', end_str='', total_length=0):
#     bar = ''.join(['#'] * int(percent * total_length)) + ''
#     bar = '\r' + start_str + bar.ljust(total_length) + ' {:0>4.1f}%|'.format(percent*100) + end_str
#     print(bar, end='', flush=True)
#
# def input_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     #resultFolder = "./result/"
#     resultFolder = program_folder+ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#     unique_inputs = get_unique(valid_input)
#     # print(unique_inputs)
#     #input_num = len(unique_inputs)
#     r = int(K/M)
#     counts = np.zeros((len(valid_input),M))
#     fre = np.zeros((len(valid_input),M))
#     co = 0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_input_coverage_'+module_name+'.txt','w')
#     result_file = open(resultFolder + 'RESULTS_input_coverage_'+module_name+'.txt','w')
#     ass_file = open(resultFolder + 'ASSESSMENT_input_coverage_'+module_name+'.txt','w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co/K, start_str='',end_str='100%',total_length=30)
#             for i in unique_inputs: #i为input
#                 count_cases += 1
#                 start = time.time()
#                 result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                 end = time.time()
#                 T2 += end-start
#                 input_file.write(str(i)+' ')
#                 # print(i)
#                 # print(result)
#                 o = list(result)[0]
#                 result_file.write('('+str(i)+','+str(o)+')'+' ')
#                 if check_WOO(i, o, valid_input, valid_output) == False:
#                     ass_file.write("Fail for WOO")
#                     result_file.write('\n')
#                     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                     result_file.write('\n')
#                     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                     # print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                     # print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                     return
#                     #end_running()
#                 for mark in range(len(valid_input)):
#                     if valid_input[mark] == i and valid_output[mark] == o:
#                         counts[mark][g] += 1
#             input_file.write('\n')
#             result_file.write('\n')
#
#     result_file.write('The total number of test cases is '+ str(count_cases)+'.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     #getting frequency
#     sum =np.zeros(M)
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i+1]
#         # print("start=" + str(start))
#         # print("end=" + str(end))
#         for j in range(start, end):
#             sum += counts[j]
#         for j in range(start, end):
#             fre[j] = counts[j]/sum
#         # print(sum)
#         sum = np.zeros(M)
#
#     #print(fre)
#     pvalue = wilcoxon(fre,p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
#
# def input_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     #input_num = len(unique_inputs)
#     r = int(K/M)
#     counts = np.zeros((len(valid_input),M))
#     fre = np.zeros((len(valid_input),M))
#     co = 0
#     count_cases = 0
#     sum = np.zeros((len(all_inputs),M))
#     input_full, output_full = check_partial_ps(valid_input, valid_output, p)
#     # print("input full:"+str(input_full))
#     # print("output full:"+str(output_full))
#
#     input_file = open(resultFolder + 'INPUTS_input_coverage_'+module_name+'.txt','w')
#     result_file = open(resultFolder + 'RESULTS_input_coverage_'+module_name+'.txt','w')
#     ass_file = open(resultFolder + 'ASSESSMENT_input_coverage_'+module_name+'.txt','w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             for i in all_inputs: #i为input
#                 count_cases += 1
#                 start = time.time()
#                 result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                 end = time.time()
#                 T2 += end - start
#                 input_file.write(str(i) + ' ')
#                 sum[int(i,2)][g] += 1
#                 # print(i)
#                 # print(result)
#                 o = list(result)[0]
#                 result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                 if i in input_full:
#                     if check_WOO(i, o, input_full, output_full) == False:
#                         ass_file.write("Fail for WOO")
#                         result_file.write('\n')
#                         result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                         result_file.write('\n')
#                         result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                         print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         end_running()
#                 #print('input='+i)
#                 #print('output='+o)
#                 for mark in range(len(valid_input)):
#                     if valid_input[mark] == i and valid_output[mark] == o:
#                         counts[mark][g] += 1
#             input_file.write('\n')
#             result_file.write('\n')
#
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2)+"s.")
#
#     #getting frequency
#     #sum =np.zeros(M)
#     # print(counts)
#     # print(sum)
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i+1]
#         for j in range(start, end):
#             fre[j] = counts[j]/sum[int(valid_input[start],2)]
#         # print(sum)
#
#     # print(fre)
#     pvalue = wilcoxon(fre, p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
#
# def input_coverage_no(inputID, outputID, num_qubit, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     r = int(K/M)
#     co = 0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_input_coverage_'+module_name+'.txt','w')
#     result_file = open(resultFolder + 'RESULTS_input_coverage_'+module_name+'.txt','w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             for i in all_inputs: #i为input
#                 count_cases += 1
#                 start = time.time()
#                 result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                 end = time.time()
#                 T2 += end - start
#                 input_file.write(str(i) + ' ')
#                 # print(i)
#                 # print(result)
#                 o = list(result)[0]
#                 result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#             input_file.write('\n')
#             result_file.write('\n')
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     input_file.close()
#     result_file.close()
#
#
#
# def output_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     unique_outputs = get_unique(valid_output)
#     unique_inputs = get_unique(valid_input)
#     # print("valid_input:"+str(valid_input))
#     # print("unique_input:"+str(unique_inputs))
#     flag = np.zeros(len(unique_outputs))
#     r = int(K/M)
#     counts = np.zeros((len(valid_input),M))
#     fre = np.zeros((len(valid_input),M))
#     stop = False
#     co=0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_output_coverage_'+module_name+'.txt','w')
#     result_file = open(resultFolder + 'RESULTS_output_coverage_'+module_name+'.txt','w')
#     ass_file = open(resultFolder + 'ASSESSMENT_output_coverage_'+module_name+'.txt','w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             exe_count = 0
#             while True:
#                 for i in unique_inputs:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end - start
#                     input_file.write(str(i) + ' ')
#                     # print(result)
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     exe_count += 1
#                     if check_WOO(i, o , valid_input, valid_output) == False:
#                         ass_file.write("Fail for WOO")
#                         result_file.write('\n')
#                         result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                         result_file.write('\n')
#                         result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                         print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         end_running()
#                     for mark in range(len(valid_input)):
#                         if valid_input[mark] == i and valid_output[mark] == o:
#                             counts[mark][g] += 1
#                     if flag[unique_outputs.index(o)] == 0:
#                         flag[unique_outputs.index(o)] = 1
#                     if check_same(flag,1):
#                         stop = True
#                         break
#                     if exe_count == BUDGET:
#                         stop = True
#                         break
#                 if stop == True:
#                     flag = np.zeros(len(unique_outputs))
#                     stop = False
#                     break
#             input_file.write('\n')
#             result_file.write('\n')
#
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     #getting frequency
#     sum =np.zeros(M)
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i+1]
#         # print("start=" + str(start))
#         # print("end=" + str(end))
#         for j in range(start, end):
#             sum += counts[j]
#         for j in range(start, end):
#             fre[j] = counts[j]/sum
#         # print(sum)
#         sum = np.zeros(M)
#
#     # print(fre)
#     pvalue = wilcoxon(fre, p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
# def output_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     all_outputs = get_all(len(outputID))
#     r = int(K / M)
#     stop = False
#     co = 0
#     count_cases = 0
#
#     counts = np.zeros((len(valid_input),M))
#     flag = np.zeros(len(all_outputs))
#     fre = np.zeros((len(valid_input),M))
#     input_full, output_full = check_partial_ps(valid_input, valid_output, p)
#     sum = np.zeros((len(all_inputs),M))
#     # print("input full:"+str(input_full))
#     # print("output full:"+str(output_full))
#
#     input_file = open(resultFolder + 'INPUTS_output_coverage_'+module_name+'.txt', 'w')
#     result_file = open(resultFolder + 'RESULTS_output_coverage_'+module_name+'.txt', 'w')
#     ass_file = open(resultFolder + 'ASSESSMENT_output_coverage_'+module_name+'.txt', 'w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             exe_count = 0
#             while True:
#                 for i in all_inputs:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end - start
#                     input_file.write(str(i) + ' ')
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     sum[int(i,2)][g] += 1
#                     if i in input_full:
#                         if check_WOO(i, o, input_full, output_full) == False:
#                             ass_file.write("Fail for WOO")
#                             result_file.write('\n')
#                             result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                             result_file.write('\n')
#                             result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                             print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                             print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                             end_running()
#                     # print(i)
#                     # print(result)
#                     exe_count += 1
#                     for mark in range(len(valid_input)):
#                         if valid_input[mark] == i and valid_output[mark] == o:
#                             counts[mark][g] += 1
#                     if flag[all_outputs.index(o)] == 0:
#                         flag[all_outputs.index(o)] = 1
#                     if check_same(flag,1):
#                         stop = True
#                         break
#                     if exe_count == BUDGET:
#                         stop = True
#                         break
#                 if stop == True:
#                     stop = False
#                     flag = np.zeros(len(all_outputs))
#                     break
#             input_file.write('\n')
#             result_file.write('\n')
#
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     # print(sum)
#     #getting frequency
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i+1]
#         # print("start=" + str(start))
#         # print("end=" + str(end))
#         for j in range(start, end):
#             fre[j] = counts[j]/sum[int(valid_input[start],2)]
#         # print(sum)
#
#     # print(fre)
#     pvalue = wilcoxon(fre, p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
# def output_coverage_no(inputID, outputID, num_qubit, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     all_outputs = get_all(len(outputID))
#     flag = np.zeros(len(all_outputs))
#     r = int(K/M)
#     stop = False
#     co=0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_output_coverage_'+module_name+'.txt', 'w')
#     result_file = open(resultFolder + 'RESULTS_output_coverage_'+module_name+'.txt', 'w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             exe_count = 0
#             while True:
#                 for i in all_inputs:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end - start
#                     input_file.write(str(i) + ' ')
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     # print(i)
#                     # print(result)
#                     exe_count += 1
#                     if flag[all_outputs.index(o)] == 0:
#                         flag[all_outputs.index(o)] = 1
#                     if check_same(flag,1):
#                         stop = True
#                         break
#                     if exe_count == BUDGET:
#                         stop = True
#                         break
#                 if stop == True:
#                     stop = False
#                     flag = np.zeros(len(all_outputs))
#                     break
#             input_file.write('\n')
#             result_file.write('\n')
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#
#
# def input_output_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     unique_inputs = get_unique(valid_input)
#     BUDGET_I = int(BUDGET / len(unique_inputs))
#     flag = np.zeros(len(valid_input))
#     r = int(K / M)
#     counts = np.zeros((len(valid_input), M))
#     fre = np.zeros((len(valid_input), M))
#     input_index = input_group(valid_input)
#     co = 0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_input_output_coverage_'+module_name+'.txt', 'w')
#     result_file = open(resultFolder + 'RESULTS_input_output_coverage_'+module_name+'.txt', 'w')
#     ass_file = open(resultFolder + 'ASSESSMENT_input_output_coverage_'+module_name+'.txt', 'w')
#
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             i_index = 0
#             for i in unique_inputs:
#                 exe_count = 0 #count of executing of one input
#                 while True:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end - start
#                     input_file.write(str(i) + ' ')
#                     # print(i)
#                     # print(result)
#                     exe_count += 1
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     if check_WOO(i, o , valid_input, valid_output) == False:
#                         ass_file.write("Fail for WOO")
#                         result_file.write('\n')
#                         result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                         result_file.write('\n')
#                         result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                         print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                         end_running()
#                     for mark in range(len(valid_input)):
#                         if valid_input[mark] == i and valid_output[mark] == o:
#                             counts[mark][g] += 1
#                             break
#                     if flag[mark] == 0:
#                         flag[mark] = 1
#                     start = input_index[i_index]
#                     if i_index == len(input_index) - 1:
#                         end = len(valid_input)
#                     else:
#                         end = input_index[i_index+1]
#                     if check_same(flag[start:end], 1) == True:
#                         break
#                     if exe_count == BUDGET_I:
#                         break
#                 i_index += 1
#             input_file.write('\n')
#             result_file.write('\n')
#             flag = np.zeros(len(valid_input))
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     # getting frequency
#     sum = np.zeros(M)
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i + 1]
#         # print("start=" + str(start))
#         # print("end=" + str(end))
#         for j in range(start, end):
#             sum += counts[j]
#         for j in range(start, end):
#             fre[j] = counts[j] / sum
#         # print(sum)
#         sum = np.zeros(M)
#
#     # print(fre)
#     pvalue = wilcoxon(fre, p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
# def input_output_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     all_outputs = get_all(len(outputID))
#     BUDGET_I = int(BUDGET / len(all_inputs))
#     counts = np.zeros((len(valid_input), M))
#     fre = np.zeros((len(valid_input), M))
#     flag = np.zeros(len(all_outputs))
#     r = int(K / M)
#     co = 0
#     count_cases = 0
#     input_full, output_full = check_partial_ps(valid_input, valid_output, p)
#     sum = np.zeros((len(all_inputs), M))
#     # print("input full:"+str(input_full))
#     # print("output full:"+str(output_full))
#
#     input_file = open(resultFolder + 'INPUTS_input_output_coverage_'+module_name+'.txt', 'w')
#     result_file = open(resultFolder + 'RESULTS_input_output_coverage_'+module_name+'.txt', 'w')
#     ass_file = open(resultFolder + 'ASSESSMENT_input_output_coverage_'+module_name+'.txt', 'w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             for i in all_inputs:
#                 exe_count = 0 #count of executing of one input
#                 while True:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end-start
#                     input_file.write(str(i) + ' ')
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     sum[int(i, 2)][g] += 1
#                     exe_count += 1
#                     if i in input_full:
#                         if check_WOO(i, o, input_full, output_full) == False:
#                             ass_file.write("Fail for WOO")
#                             result_file.write('\n')
#                             result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#                             result_file.write('\n')
#                             result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                             print("The total run time is " + "{:.2f}".format(T1) + "s.")
#                             print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#                             end_running()
#                     if flag[all_outputs.index(o)] == 0:
#                         flag[all_outputs.index(o)] = 1
#                     for mark in range(len(valid_input)):
#                         if valid_input[mark] == i and valid_output[mark] == o:
#                             counts[mark][g] += 1
#                     if check_same(flag,1):
#                         break
#                     if exe_count == BUDGET_I:
#                         break
#                 flag = np.zeros(len(all_outputs))
#             input_file.write('\n')
#             result_file.write('\n')
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#     # print(sum)
#     # getting frequency
#     input_index = input_group(valid_input)
#     # print(input_index)
#     for i in range(len(input_index)):
#         start = input_index[i]
#         if i == len(input_index) - 1:
#             end = len(valid_input)
#         else:
#             end = input_index[i + 1]
#         # print("start=" + str(start))
#         # print("end=" + str(end))
#         for j in range(start, end):
#             fre[j] = counts[j] / sum[int(valid_input[start], 2)]
#         # print(sum)
#
#     # print(fre)
#     pvalue = wilcoxon(fre, p)
#     judge_ass_result(valid_input, valid_output, pvalue, ass_file)
#
#     input_file.close()
#     result_file.close()
#     ass_file.close()
#
# def input_output_coverage_no(inputID, outputID, num_qubit, module_name, program_folder):
#     global T2
#     resultFolder = program_folder + ROOT+"result"+ROOT
#     if not os.path.exists(resultFolder):
#         os.makedirs(resultFolder)
#
#     all_inputs = get_all(len(inputID))
#     all_outputs = get_all(len(outputID))
#     BUDGET_I = int(BUDGET / len(all_inputs))
#     flag = np.zeros(len(all_outputs))
#     r = int(K / M)
#     co = 0
#     count_cases = 0
#
#     input_file = open(resultFolder + 'INPUTS_input_output_coverage_'+module_name+'.txt', 'w')
#     result_file = open(resultFolder + 'RESULTS_input_output_coverage_'+module_name+'.txt', 'w')
#
#     for g in range(M):
#         for t in range(r):
#             co += 1
#             process_bar(co / K, start_str='', end_str='100%', total_length=30)
#             for i in all_inputs:
#                 exe_count = 0 #count of executing of one input
#                 while True:
#                     count_cases += 1
#                     start = time.time()
#                     result = execute_quantum_program(inputID, outputID, num_qubit, i, module_name)
#                     end = time.time()
#                     T2 += end-start
#                     input_file.write(str(i) + ' ')
#                     o = list(result)[0]
#                     result_file.write('(' + str(i) + ',' + str(o) + ')' + ' ')
#                     # print(i)
#                     # print(result)
#                     exe_count += 1
#                     if flag[all_outputs.index(o)] == 0:
#                         flag[all_outputs.index(o)] = 1
#                     if check_same(flag, 1) == True:
#                         break
#                     if exe_count == BUDGET_I:
#                         break
#                 flag = np.zeros(len(all_outputs))
#             input_file.write('\n')
#             result_file.write('\n')
#     result_file.write('The total number of test cases is ' + str(count_cases) + '.')
#     result_file.write('\n')
#     result_file.write('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")
#
# def check_configuration_file(config):
#     if config.has_section('program') == False:
#         print("Error: Quito cannot find section 'program' in this configuration file.")
#         end_running()
#     else:
#         if config.has_option('program', 'root') == False:
#             print("Error: Quito cannot find the root of the program.")
#             end_running()
#         if config.has_option('program', 'num_qubit') == False:
#             print("Error: Quito cannot find the number of qubits of the program.")
#             end_running()
#         if config.has_option('program', 'inputID') == False:
#             print("Error: Quito cannot find the input IDs of the program.")
#             end_running()
#         if config.has_option('program', 'outputID') == False:
#             print("Error: Quito cannot find the output IDs of the program.")
#             end_running()
#     if config.has_section('program_specification_category') == False:
#         print("Error: Quito cannot find section 'program_specification_category' in this configuration file.")
#         end_running()
#     else:
#         if config.has_option('program_specification_category', 'ps_category') == False:
#             print("Error: Quito cannot find the category of the program specification.")
#             end_running()
#     if config.has_section('quito_configuration') == False:
#         print("Error: Quito cannot find section 'quito_configuration' in this configuration file.")
#         end_running()
#     else:
#         if config.has_option('quito_configuration', 'coverage_criterion') == False:
#             print("Error: Quito cannot find the coverage criterion you choose.")
#             end_running()
#
#     ps_category = config.get('program_specification_category', 'ps_category')
#     if ps_category == 'full' or ps_category == 'partial':
#         if config.has_section('program_specification') == False:
#             print("Error: Quito cannot find the program specification.")
#             end_running()
#     return ps_category
#
# def check_inputID_outputID(num_qubit, inputID, outputID):
#     if check_unique(inputID) == False:
#         raise Exception("Wrong input IDs")
#     if check_unique(outputID) == False:
#         raise Exception("Wrong output IDs")
#     inputID.sort()
#     outputID.sort()
#
#     if int(inputID[-1]) > num_qubit - 1:
#         raise Exception("Wrong input IDs")
#     if int(inputID[-1]) > num_qubit - 1:
#         raise Exception("Wrong output IDs")
#
#     return inputID, outputID
#
# def check_bin(bin_str, n):
#     if len(bin_str) != n:
#         print("Error: The format of the program specification is wrong.")
#         end_running()
#    # print("check bin: "+str(bin_str))
#     for i in range(len(bin_str)):
#         if bin_str[i] != '0' and bin_str[i] != '1':
#             print("Error: The format of the program specification is wrong.")
#             end_running()
#
# def check_parameters(inputID, outputID, num_qubit, module_name, program_folder):
#     if isinstance(inputID, (list, numpy.ndarray)) == False:
#         raise Exception(f'Invalid parameter type {type(inputID)} for inputID')
#     if isinstance(outputID, (list, numpy.ndarray)) == False:
#         raise Exception(f'Invalid parameter type {type(outputID)} for outputID')
#     if isinstance(num_qubit, int) == False:
#         raise Exception(f'Invalid parameter type {type(num_qubit)} for num_qubit')
#     inputID, outputID = check_inputID_outputID(num_qubit, inputID, outputID)
#
#
#
# def quito_run():
#     global START
#     print("please enter the root of your configuration file.(.ini file)")
#
#     #get configuration file
#     root_con = input()
#     START = time.time()
#     if os.path.isfile(root_con) == True:
#         config = configparser.ConfigParser(allow_no_value=True)
#         config.read(root_con, encoding='utf-8')
#     else:
#         print("Error: Quito cannot find the configuration file.")
#         end_running()
#
#     ps_category = check_configuration_file(config)
#     if ps_category != 'no' and ps_category != 'partial' and ps_category != 'full':
#         print("Error: The format of program specification category is wrong.")
#         end_running()
#     print(ps_category)
#     #get quantum program
#     root = config.get('program','root')
#     # print(root)
#
#     if os.path.isfile(root) != True:
#         print("Error: Quito cannot find the quantum program file.")
#         end_running()
#
#     root_list = root.split(ROOT)
#     program_file = root_list[len(root_list)-1]
#     program_folder = root_list[:len(root_list)-1]
#     program_folder = ROOT.join(str(i) for i in program_folder)
#     sys.path.append(program_folder)
#     # print(program_file.split('.')[0])
#     module_name = program_file.split('.')[0]
#
#     #get inputID, outputID and numner of qubits
#     inputID_o = config.get('program','inputID').split(',')
#     outputID_o = config.get('program','outputID').split(',')
#     num_qubit = int(config.get('program','num_qubit'))
#     inputID, outputID = check_inputID_outputID(num_qubit, inputID_o, outputID_o)
#
#     #ps_category = config.get('program_specification_category','ps_category')
#     coverage_criterion = config.get('quito_configuration', 'coverage_criterion')
#     print(coverage_criterion)
#     if coverage_criterion != 'IC' and coverage_criterion != 'OC' and coverage_criterion != 'IOC':
#         print("Error: The format of coverage criterion is not right.")
#         end_running()
#
#     if config.get('quito_configuration', 'K') != None:
#         global K
#         K = int(config.get('quito_configuration', 'K'))
#     if config.get('quito_configuration', 'M') != None:
#         global  M
#         M = int(config.get('quito_configuration', 'M'))
#     if config.get('quito_configuration', 'confidence_level') != None:
#         global C_LEVEL
#         C_LEVEL = float(config.get('quito_configuration', 'confidence_level'))
#     if config.get('quito_configuration', 'statistical_test') != None:
#         global  S_TEST
#         S_TEST = config.get('quito_configuration', 'statistical_test')
#
#     global BUDGET
#
#
#     # print('num_qubit:'+str(num_qubit))
#
#
#
#     if ps_category == 'no':
#         if coverage_criterion == 'IC':
#             input_coverage_no(inputID, outputID, num_qubit, module_name, program_folder)
#         else:
#             BUDGET = pow(2, len(inputID)) * 10  # default
#             if config.get('quito_configuration', 'BUDGET') != None:
#                 BUDGET = int(config.get('quito_configuration', 'BUDGET'))
#
#             # check budget
#             if BUDGET < pow(2, len(inputID)):
#                 print("Error: Budget is smaller than the number of inputs.")
#                 end_running()
#             if coverage_criterion == 'OC':
#                 output_coverage_no(inputID, outputID, num_qubit, module_name, program_folder)
#             elif coverage_criterion == 'IOC':
#                 input_output_coverage_no(inputID, outputID, num_qubit, module_name, program_folder)
#
#
#     else:
#         #get PS
#         valid_input = []
#         valid_output = []
#         p = []
#         ps = config.items('program_specification')
#         #print("origin:"+str(ps))
#
#         if check_unique(ps) == False:
#             print("Program specifications not unique")
#             end_running()
#         #sort PS according to input and output
#         ps.sort(key=lambda x:x[0])
#         #print("new:"+str(ps))
#         for i in range(len(ps)):
#             valid_input_item = ps[i][0][:len(inputID)]
#             valid_output_item = ps[i][0][len(inputID)+1:]
#             check_bin(valid_input_item,len(inputID))
#             check_bin(valid_output_item,len(outputID))
#             valid_input.append(valid_input_item)
#             valid_output.append(valid_output_item)
#             p.append(float(ps[i][1]))
#
#         # print(valid_input)
#         # print(valid_output)
#         # print(p)
#
#
#         # check budget
#         if coverage_criterion == 'OC' or 'IOC':
#             if ps_category == 'full':
#                 BUDGET = pow(2, len(set(valid_input))) * 10  # default
#                 if config.get('quito_configuration', 'BUDGET') != None:
#                     BUDGET = int(config.get('quito_configuration', 'BUDGET'))
#                 num_valid_input = set(valid_input)
#                 if BUDGET < len(num_valid_input):
#                     print("Error: Budget is smaller than the number of inputs.")
#                     end_running()
#             elif ps_category == 'partial':
#                 BUDGET = pow(2, len(inputID)) * 10  # default
#                 if config.get('quito_configuration', 'BUDGET') != None:
#                     BUDGET = int(config.get('quito_configuration', 'BUDGET'))
#                 if BUDGET < pow(2,len(inputID)):
#                     print("Error: Budget is smaller than the number of inputs.")
#                     end_running()
#
#         if ps_category == 'full':
#             check_full_ps(valid_input, p)
#         # print(outputID)
#
#         if ps_category == 'full':
#             if coverage_criterion == 'IC':
#                 input_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#             elif coverage_criterion == 'OC':
#                 output_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#             elif coverage_criterion == 'IOC':
#                 input_output_coverage(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#
#         if ps_category == 'partial':
#             if coverage_criterion == 'IC':
#                 input_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#             elif coverage_criterion == 'OC':
#                 output_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#             elif coverage_criterion == 'IOC':
#                 input_output_coverage_partial(inputID, valid_input, valid_output, num_qubit, outputID, p, module_name, program_folder)
#
# def user_input():
#     print('\n')
#     text1 = "1. Check the template of the configuration file.(.ini file)"
#     print(text1)
#     print("\n")
#     text2 = "2. Check the example of the configuration file.(.ini file)"
#     print(text2)
#     print('\n')
#     text3 = "3. Upload your configuration.(.ini file)"
#     print(text3)
#     print('\n')
#     print("Your choice?[1/2/3]",end=' ')
#     user_choice = input()
#     print('\n')
#     return user_choice
#
# def template_file():
#     print("The template is shown below:")
#     program = '''
# [program]
# root=
# ;(Required)
# ;Description: The absolute root of your quantum program file.
# num_qubit=
# ;(Required)
# ;Description: The total number of qubit of your quantum program.
# inputID=
# ;(Required)
# ;Description: The ID of input qubits.
# outputID=
# ;(Required)
# ;Description: The ID of output qubits which are the qubits to be measured.
# ;Format: A non-repeating sequence separated by commas.
#     '''
#     quito_configuration = '''
# [quito_configuration]
# coverage_criterion=
# ;Description: The coverage criterion you choose.
# ;Choice: IC/OC/IOC
# K=
# ;(Optional)
# ;Description: The total number of test suites, K=200 by default.
# M=
# ;(Optional)
# ;Description: The number of test suite groups, M=20 by default.
# BUDGET=
# ;(Optional)
# ;Description: The budget of the number of test cases in one test suite, BUDGET=10*number of inputs by default.
# confidence_level=
# ;(Optional)
# ;Description: The confidence level for statistical test, confidence_level=0.01 by default.
# statistical_test=
# ;(Optional)
# ;Description: The statistical test for assessment, statistical_test=one-sample Wilcoxon signed rank test by default.
# '''
#     program_specification_category = '''
# [program_specification_category]
# ps_category=
# ;(Required) Description: The category of your program specification. (full/partial/no)
# '''
#     program_specification = '''
# [program_specification]
# ;(Required for full and partial program specification)
# ;Description: The program specification.
# ;Format:input string,output string=probability
# '''
#     print(program)
#     print(quito_configuration)
#     print(program_specification_category)
#     print(program_specification)
#     print('\n')
#
# def example_file():
#     print("The example is shown below:")
#     example = '''
# [program]
# root=C:\program_run.py
# num_qubit=3
# inputID=0,1
# outputID=2
#
# [program_specification_category]
# ps_category=full
#
# [quito_configuration]
# coverage_criterion=IC
# K=200
# M=10
# BUDGET=20
# confidence_level=0.01
# statistical_test=one-sample Wilcoxon signed rank test
#
# [program_specification]
# 01,1=0.5
# 01,0=0.5
# 00,1=1
# 11,1=1
# 10,1=0.5
# 10,0=0.5
# '''
#     print(example)
#
# if __name__ == '__main__':
#     warnings.filterwarnings('ignore')
#     print("".center(60,"="))
#     print('\n')
#     print("Welcome to Quito".center(60))
#     print('\n')
#     print("".center(60,"="))
#     user_choice = user_input()
#     while user_choice != '3':
#         if user_choice != '1' and user_choice != '2':
#             print("Error: Please type in '1', '2' or '3'.")
#             user_choice = user_input()
#         elif user_choice == '1':
#             print("template")
#             template_file()
#             user_choice = user_input()
#         elif user_choice == '2':
#             print("example")
#             example_file()
#             user_choice = user_input()
#     if user_choice == '3':
#         #start = time.time()
#         quito_run()
#         end = time.time()
#         T1 = end - START
#         print('\n')
#         print("The total run time is "+"{:.2f}".format(T1)+"s.")
#         print('The execution time of the quantum program is ' + "{:.2f}".format(T2) + "s.")

from quito import quito

quito("/home/xinyi/quito/ex3/quito.py")

