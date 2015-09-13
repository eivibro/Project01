#Python file for plotting results
import os 

from pylab import *
with open('exact_and_approximated_n100.txt', 'r') as f:
	line = f.readline()
	data0 = line.split()
	number_of_lines = int(data0[0])
	exact = zeros(number_of_lines)
	estimated = zeros(number_of_lines)
	error = zeros(number_of_lines)
	for i in range(0, number_of_lines - 1):
		line = f.readline()
		data = line.split()	
		exact[i] = data[0]
		estimated[i] = data[1]
		error[i] = data[2]

x = zeros(number_of_lines)
h = 1./(number_of_lines-1)
for i in range(0,number_of_lines):
	x[i] = i*h
print max(error[2:len(error)-1])
plot(x, exact, '.-', x, estimated, '.-')
xlabel('$x$', fontsize = 20)
ylabel('$u(x)$', fontsize = 20)
xticks(fontsize = 14) 
yticks(fontsize = 14)
legend(['Exact', 'Approximated'], fontsize = 16)
filename = 'second_derivative_n'+ str(number_of_lines)+ '.eps'
print filename
savefig(filename)
cmd = "mv "+filename + " ../latex/."
os.system(cmd)
os.system("cp ../Project01/main.cpp ../latex/.")
os.system("cp ../build-Project01-Desktop-Debug/reader.py ../latex/.")
show()

#Plotting the errors as a function of the step length
with open('errors_for_different_steplengths.txt', 'r') as f:
	line = f.readline()
	data0 = line.split()
	number_of_lines = int(data0[0])
	steplength = zeros(number_of_lines)
	errors_for_different_steplengths = zeros(number_of_lines)
	for i in range(0, number_of_lines):
		line = f.readline()
		data = line.split()
		steplength[i] = data[0]
		errors_for_different_steplengths[i] = data[1]

loglog(steplength, 10**errors_for_different_steplengths)
xlabel('step length', fontsize = 14)
ylabel('Error', fontsize = 14)
xticks(fontsize = 12) 
yticks(fontsize = 12)
savefig('Error_func.eps')
os.system('cp Error_func.eps ../latex')
show()