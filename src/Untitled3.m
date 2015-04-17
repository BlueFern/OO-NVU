
prompt = {'Enter size of population:',
            'Enter number of organism which are allowed to reproduce:',
            'Enter lower boundary of randomisation of the initial population:',
            'Enter upper boundary of randomisation of the initial population:',
            'Enter the probability of mutation:',
            'Enter lower boundary of mutation:',
            'Enter upper boundary of mutation:',
            'Enter the probability of cross-over:',
            'Enter the threshold for elitism:'};
dlg_title = 'Parameter sensitivity test';
num_lines = 1;
defaultanswer = {'30','10','0.8','1.2','0.1','0.8','1.2','0.8','0.005'};

x = inputdlg(prompt,dlg_title,num_lines,defaultanswer);

num_pop     = str2num(x{1});
nr_parents  = str2num(x{2});
low_bound 	= str2num(x{3});
up_bound    = str2num(x{4});
p_mu        = str2num(x{5});
mu_low_b    = str2num(x{6});
mu_up_b     = str2num(x{7});
p_co        = str2num(x{8});
p_el        = str2num(x{9});

 
