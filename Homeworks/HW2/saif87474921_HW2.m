%%
% <latex>
% \title{BE 521: Homework 2 Questions \\{\normalsize Modeling Neurons} \\{\normalsize Spring 2021}}
% \author{46 points}
% \date{Due: Tuesday, 2/2/2021 10:00 PM}
% \maketitle
% \textbf{Objective:} Computational modeling of neurons. \\
% We gratefully acknowledge Dr. Vijay Balasubramanian (UPenn) for many of
% the questions in this homework.\\
% </latex>

%% 
% <latex>
% \begin{center}
% \author{Saif Khawaja \\
%   \normalsize Collaborators: Raveen K \\}
% \end{center}
% </latex>

%%
% <latex>
% \section{Basic Membrane and Equilibrium Potentials (6 pts)}
% Before undertaking this section, you may find it useful to read pg.
% 153-161 of Dayan \& Abbott's \textit{Theoretical Neuroscience} (the 
% relevant section of which, Chapter 5, is posted with the homework). 
% </latex>

%%
% <latex>
% \begin{enumerate}
%  \item Recall that the potential difference $V_T$ when a mole of ions crosses a cell membrane is defined by the universal gas constant $R = 8.31\; {\rm J/mol\cdot K}$, the temperature $T$ (in Kelvin), and Faraday's constant $F = 96,480 {\rm\ C/mol}$ \[ V_T = \frac{RT}{F} \] Calculate $V_T$ at human physiologic temperature ($37\; ^{\circ} \rm C$). (1 pt)
% </latex>

%%
% <latex>
% \rule{\textwidth}{1pt}
% \textit{Example Latex math commands that uses the align tag to make your equations
% neater. You can also input math into sentences with \$ symbol like $\pi + 1$.}
% \begin{align*}
% E = MC^2 \tag{not aligned}\\
% E = & MC^2 \tag{aligned at = by the \&}\\
% 1 = &\; \frac{2}{2}\tag{aligned at = by \&}
% \end{align*}
% \rule{\textwidth}{1pt}
% </latex>

Vt = (8.31*310.15)/96480 *1000

% Vt at the temperature of the human body is 26.7138 .
%%
% <latex>
%  \item Use this value $V_T$ to calculate the Nernst equilibrium potentials 
%  (in mV) for the $\rm K^+$, $\rm Na^+$, and $\rm Cl^-$ ions, given the following 
%  cytoplasm and extracellular concentrations in the squid giant axon: 
%  $\rm K^+$ : (120, 4.5), $\rm Na^+$ : (15, 145), and $\rm Cl^-$ : (12, 120), 
%  where the first number is the cytoplasmic and the second the extracellular 
%  concentration (in mM). (2 pts)
% </latex>

conc_K = [120, 4.5];
conc_Na = [15, 145];
conc_Cl = [12, 120];

R = 8.31; 
T = 310.15; 
F = 96480;

Nernst_K = ((R * T) / (1 * F) * log(conc_K(2) / conc_K(1))) * 1000
Nernst_Na = ((R * T) / (1 * F) * log(conc_Na(2) / conc_Na(1))) * 1000
Nernst_Cl = ((R * T) / -(1 * F) * log(conc_Cl(2) / conc_Cl(1))) * 1000

% For K, Na, and Cl respectively, the Nernst potentials are -87.7124,
% 60.6051, and -61.5108 mV.

%%
% <latex>
%  \item 
%   \begin{enumerate}
% 	\item Use the Goldmann equation,
% 	  \begin{equation}
% 		V_m = V_T\ln\left( \frac{\rm P_{K}\cdot[K^+]_{out} + P_{NA}\cdot[Na^+]_{out} + P_{Cl}\cdot[Cl^-]_{in}}{\rm P_{K}\cdot[K^+]_{in} + P_{NA}\cdot[Na^+]_{in} + P_{Cl}\cdot[Cl^-]_{out}} \right)
% 	  \end{equation}
% 	to calculate the resting membrane potential, $V_m$, assuming that the ratio of membrane permeabilities $\rm P_K:P_{Na}:P_{Cl}$ is $1.0:0.045:0.45$. Use the ion concentrations given above in Question 1.2. (2 pts)
% </latex>

permeabilities = [1, 0.05, 0.45];

membrane_potential = Vt * log((permeabilities(1) * conc_K(2) + permeabilities(2) * conc_Na(2) + permeabilities(3) * conc_Cl(1))/(permeabilities(1) * conc_K(1) + permeabilities(2) * conc_Na(1) + permeabilities(3) * conc_Cl(2)))

% The resting membrane potential is -62.0123 mV.

%%
% <latex>
% 	\item Calculate the membrane potential at the peak action potential, assuming a permeability ratio of $1.0:11:0.45$, again using the ion concentrations given in Question 1.2. (1 pt)
%   \end{enumerate}
% </latex>

permeabilities = [1, 12, 0.45];

membrane_potential_peak = Vt * log((permeabilities(1) * conc_K(2) + permeabilities(2) * conc_Na(2) + permeabilities(3) * conc_Cl(1)) / (permeabilities(1) * conc_K(1) + permeabilities(2) * conc_Na(1) + permeabilities(3) * conc_Cl(2)))

% The membrane potential at peak action potential is 42.6891 mV.

%%
% <latex>
% 	\item The amplitudes of the multi-unit signals in HW0 and local field
% 	potentials (LFP) in HW1 had magnitudes on the order of 10 to 100
% 	microvolts. The voltage at the peak of the action potential (determined
% 	using the Goldman equation above) has a magnitude on the order of 10
% 	millivolts. Briefly explain why we see this difference in magnitude.
% 	Hint 1: Voltage is the difference in electric potential between two
% 	points. What are the two points for our voltage measurement in the
% 	multi-unit and LFP signals? What are the two points for the voltage
% 	measurement of the action potential? Hint 2: The resistance of the neuronal membrane is typically much higher than the resistance of the extracellular fluid. (2 pts)
% </latex>

% Since voltage is proportional to resistance (V = IR), and the
% extracellular membrane resistance is typically much lower than the
% neuronal membrane resistance (hint 2), the peak action potential voltage
% is much higher than any of the LFP or multi-unit activity signals from
% previous homework. Furthermore, LFP readings were measured in
% extra-cellular fluid but action potentials were measured across the
% neuronal cell and, given hint 1, this influences the recorded reading. 

%% 
% <latex>
% \end{enumerate}
% \section{Integrate and Fire Model (38 pts)}
% You may find it useful to read pg.\ 162-166 of Dayan and Abbott for this section. The general differential equation for the integrate and fire model is
% \[ \tau_m\frac{dV}{dt} = V_m - V(t) + R_m I_e(t) \]
% where $\tau_m = 10\, \rm ms$ is the membrane time constant, describing how fast the current is leaking through the membrane, $V_m$ in this case is constant and represents the resting membrane potential (which you have already calculated in question 1.3.a), and $V(t)$ is the actual membrane potential as a function of time. $R_m = 10^7\, \Omega$ is the constant total membrane resistance, and $I_e(t)$ is the fluctuating incoming current. Here, we do not explicitly model the action potentials (that's Hodgkin-Huxley) but instead model the neuron's behavior leading up and after the action potential.
% </latex>

%%
% <latex>
% Use a $\Delta t = 10\, \rm \mu s$ ($\Delta t$ is the discrete analog of the continuous $dt$). Remember, one strategy for modeling differential equations like this is to start with an initial condition (here, $V(0)=V_m$), then calculate the function change (here, $\Delta V$, the discrete analog to $dV$) and then add it to the function (here, $V(t)$) to get the next value at $t+\Delta t$. Once/if the membrane potential reaches a certain threshold ($V_{th} = -50\, \rm mV$), you will say that an action potential has occurred and reset the potential back to its resting value.
% \begin{enumerate}
%  \item Model the membrane potential with a constant current injection (i.e., $I_e(t) = I_e = 2 {\rm nA}$). Plot your membrane potential as a function of time to show at least a handful of ``firings.'' (8 pts)
% </latex>

I_e = 2e-9;  

R_membrane = 1*10^7; 

V_threshold = -50e-3;

tau = 0.01;  

dt = 1e-5;  

range = 0:dt:0.1;

%%

V = [];
V(1) = (membrane_potential / 1000); 
peaks = 0;

for i=2:length(range)
    dv = (dt / tau) * ((membrane_potential / 1000) - V(i-1) + (R_membrane * I_e));
    dv;
    V(i) = dv + V(i-1);
    if V(i) >= V_threshold
        peaks = peaks + 1;
        V(i) = membrane_potential / 1000;
    end
end

V = V.*1000; 

figure
plot(range, V)
ylabel('Potential (mV)')
xlabel('T (s)')
title('Membrane Potential against Time')
ylim([-64 -49])

%%
% <latex>
%  \item Produce a plot of firing rate (in Hz) versus injection current, over the range of 1-4 nA. (4 pts)
% </latex>


I_e = 2e-9;

V_threshold=-50e-3; 

R_membrane=1*10^7; 

tau=0.01; 

dt=1e-5; 

range=0:dt:1;

V=[];
p_vector=[];
V(1)=(membrane_potential / 1000); 

currents = 1e-9:1e-10:4e-9;

for j = 1:length(currents)
    peaks = 0;
    for i = 2:length(range)
        dv = (dt/tau)*((membrane_potential / 1000) - V(i-1) + (R_membrane*currents(j)));
        dv;
        V(i) = dv+V(i-1);
        if V(i) >= V_threshold
            peaks = peaks + 1;
            V(i)= membrane_potential / 1000;
        end
        
    end
    p_vector = [p_vector peaks];
end


figure
plot(currents, p_vector)
ylabel('Rate (Hz)')
xlabel('I_e (nA)')
title('Firing Rate v. Injection Current')

%%
% <latex>
%  \item \texttt{I521\_A0002\_D001} contains a dynamic current injection in nA. Plot the membrane potential of your neuron in response to this variable injection current. Use Matlab's \texttt{subplot} function to place the plot of the membrane potential above the injection current so that they both have the same time axis. (Hint: the sampling frequency of the current injection data is different from the sampling frequency ($\frac{1}{\Delta t}$) that we used above.) (4 pts)
% </latex>

session = IEEGSession('I521_A0002_D001', 'saifkhawaja', '/Users/saif/Documents/GitHub/braincomputerinterfaces/Homeworks/sai_ieeglogin.bin' )

session.data;

durationInUSec = session.data(1).rawChannels(1).get_tsdetails.getDuration;
durationInSec = durationInUSec / 1e6;

nr = ceil((session.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session.data.sampleRate);
all_data = session.data.getvalues(1:nr, 1);
all_data = all_data.* 1e-9;

downsampled_currents = downsample(all_data, 10);
downsampled_currents = downsampled_currents.*1e-9;

sr = (1/1e-5);

I_e = 2e-9; 

V_threshold = -50e-3; 

tau = 0.01;

t = durationInSec;

dt = 1e-5; 

R_membrane = 1*10^7; 

peak_vector_i_varied = [];

t2 = 0:dt:(0.5-dt);

V_new = [];
V_new(1) = (membrane_potential / 1000);

peaks = 0;

for i=2:length(t2)
    dv=(dt / tau) * ((membrane_potential / 1000) - V_new(i-1) + (R_membrane * downsampled_currents(i)));
    dv;
    V_new(i) = dv + V_new(i-1);
    if V_new(i) >= V_threshold
        peaks = peaks+1;
        V_new(i) = membrane_potential / 1000;
    end
end

V_new = V_new.*1000;

figure
subplot(2,1,1)
plot(t2,V_new)
xlabel('Time (s)')
ylabel('Potential (mV)')
title('Membrane Potential v. Time')
lims = [-65 -48];
ylim(lims)
subplot(2,1,2)
plot(t2,downsampled_currents)
xlabel('Time (s)')
ylabel('Ie (A)')
title('Injection Current v. Time')

%%
% <latex>
%  \item Real neurons have a refractory period after an action potential that prevents them from firing again right away. We can include this behavior in the model by adding a spike-rate adaptation conductance term, $g_{sra}(t)$ (modeled as a potassium conductance), to the model
%  \[ \tau_m\frac{dV}{dt} = V_m - V(t) - r_m g_{sra}(t)(V(t)-V_K)+ R_m I_e(t) \]
%  where \[ \tau_{sra}\frac{dg_{sra}(t)}{dt} = -g_{sra}(t), \]
%  Every time an action potential occurs, we increase $g_{sra}$ by a certain constant amount, $g_{sra} = g_{sra} + \Delta g_{sra}$. Use $r_m \Delta g_{sra} = 0.06$. Use a conductance time constant of $\tau_{sra} = 100\, \rm ms$, a potassium equilibrium potential of $V_K = -70\, \rm mV$, and $g_{sra}(0) = 0$. (Hint: How can you use the $r_m \Delta g_{sra}$ value to update voltage and conductance separately in your simulation?)
% </latex>

dt = 1e-5; 

t7 = 0:dt:0.2;

I_e = 2e-9; 

tau = 0.01; 

taug = 0.1; 

R_membrane = 1*10^7; 

V_threshold = -50e-3; 

RM_increase = 0.06; 

V_k = -70e-3; 

V_4 = [];
V_4(1) = (membrane_potential/1000); 

RM_gsra = [];
RM_gsra(1) = 0;

peaks = 0;


for i=2:length(t7)
    dv_RM_gsra = (dt/taug)*(-RM_gsra(i-1));
    RM_gsra(i) = RM_gsra(i-1)+dv_RM_gsra;
    dv_4 = (dt/tau)*((membrane_potential / 1000)-V_4(i-1)-(RM_gsra(i-1) * (V_4(i-1)-V_k)) + (R_membrane*I_e));
    dv_4;
    V_4(i) = dv_4+V_4(i-1);
    if V_4(i) >= V_threshold
        RM_gsra(i) = RM_gsra(i)+RM_increase;
        V_4(i) = membrane_potential/1000;
    end
end


%%

t7 = t7.*1000;

V_4 = V_4.*1000;

%%

figure
plot(t7,V_4)
xlabel('Time (Milliseconds)')
ylabel('Membrane Potential (mV)')
ylim([-64 -49])
title('Membrane Potential v. Time')


%%
% <latex>
%  \begin{enumerate}
%   \item Implement this addition to the model (using the same other parameters as in question 2.1) and plot the membrane potential over 200 ms. (8 pts)
% </latex>

I_e = 2e-9; 

tau = 0.01; 

taug = 0.1; 

R_membrane = 1*10^7;

dt = 1e-5; 

V_threshold = -50e-3; 

RM_increase = 0.06; 

V_k = -70e-3; 

t7 = 0:dt:0.2;

RM_gsra = [];
RM_gsra(1) = 0;
V_4 = [];
V_4(1) = (membrane_potential/1000); 
peaks = 0;

for i=2:length(t7)
    dv_RM_gsra = (dt / taug) * (-RM_gsra(i-1));%R_m gsra differential eqn
    RM_gsra(i) = RM_gsra(i-1) + dv_RM_gsra;
    dv_4 = (dt/tau) * ((membrane_potential / 1000) - V_4(i-1) - (RM_gsra(i-1) * (V_4(i-1) - V_k)) + (R_membrane * I_e));
    dv_4;
    V_4(i) = dv_4+V_4(i-1);
    if V_4(i) >= V_threshold
        RM_gsra(i) = RM_gsra(i)+RM_increase;
        V_4(i) = membrane_potential/1000;
    end
end

V_4 = V_4.*1000;
t7 = t7.*1000;

figure
plot(t7, V_4)
xlabel('Time (Milliseconds)')
ylabel('Membrane Potential (mV)')
ylim([-64 -49])
title('Membrane Potential v. Time')

%%
% <latex>
%   \item Plot the inter-spike interval (the time between the spikes) of all the spikes that occur in 500 ms. (2 pts)
% </latex>

RM_gsra = [];
RM_gsra(1) = 0;
V_5 = [];
V_5(1) = (membrane_potential / 1000); 
peaks = 0;


I_e = 2e-9;

tau = 0.01; 

taug = 0.1; 

R_membrane = 1*10^7; 

dt = 1e-5;  

V_threshold = -50e-3;

RM_increase = 0.06; 

V_k = -70e-3;


t5 = 0:dt:0.5;


for i=2:length(t5)
    dv_RM_gsra = (dt/taug) * (-RM_gsra(i-1));
    RM_gsra(i) = RM_gsra(i-1) + dv_RM_gsra;
    
    dv_4=(dt/tau) * ((membrane_potential / 1000)-V_5(i - 1)-(RM_gsra(i - 1) * (V_5(i-1) - V_k)) + (R_membrane * I_e));
    dv_4;
    V_5(i) = dv_4+V_5(i-1);
    
    if V_5(i) >= V_threshold
        RM_gsra(i) = RM_gsra(i) + RM_increase;
        V_5(i) = membrane_potential/1000;
    end
end
%%

thrsh = -52;

[pks, locations] = findpeaks(V_5, 'MinPeakHeight', thrsh);

time_diff = t5(diff(locations)); %convert to seconds and find time difference

time_diff = time_diff.*1000; % convert to milliseconds

%%
figure
plot(time_diff, 'o')
ylabel('Interval (ms)')
xlabel('Count')
title('Interval v. Event Count')

%%
% <latex>
%   \item Explain how the spike-rate adaptation term we introduced above might be contributing to the behavior you observe in 2.4.b. (2 pts)
% </latex>

% Spike-rate adaption conductance supports the modelling of the period
% between potentials. This action potential refraction period uses the
% term because it increases the K conductance, which increases T between
% spiking events. This is demonstrated by the decreasing rate of the
% increase of membrane potential (dv). A greater T is needed before the
% potential can break threshold voltage and generate another spike.
% This is similar to real context where more ions create better conductance
%, which increases reset time as levels need to readjust.

%%
% <latex>
%  \end{enumerate}
%  \item Pursue an extension of this basic integrate and fire model. A few ideas are: implement the Integrate-and-Fire-or-Burst Model of Smith et al.\ 2000 (included); implement the Hodgkin-Huxley model (see Dayan and Abbot, pg.\ 173); provide some sort of interesting model of a population of neurons; or perhaps model what an electrode sampling at 200 Hz would record from the signal you produce in question 2.3. Feel free to be creative. 
%  We reserve the right to give extra credit to particularly interesting extensions and will in general be more generous with points for more difficult extensions (like the first two ideas), though it is possible to get full credit for any well-done extension.
%   \begin{enumerate}
% 	\item Briefly describe what your extension is and how you will execute it in code. (6 pts)
% </latex>

% The below simulation uses equations from the Hodgkin-Huxley model for an action potential as described by
% "Computational and Mathematical Modeling of Neural Systems" (Dayan and
% Abbott). This is a more complex model of the membrane potential when compared to the integrate and fire model. 
% Defined ICs are the foundation for computation for the for loop that
% calculates the membrane potential. From there, we can compare the
% membrane potential and probabilities of activation against time with
% plots. 

%%

permeabilities = [1, 0.05, 0.45];

membrane_potential = Vt * log((permeabilities(1) * conc_K(2) + permeabilities(2) * conc_Na(2) + permeabilities(3) * conc_Cl(1))/(permeabilities(1) * conc_K(1) + permeabilities(2) * conc_Na(1) + permeabilities(3) * conc_Cl(2)))

g = [0.003, 0.36, 1.2];
e = [-54.387, -77, 50];

cm = 0.01; 
A = 0.1; 

I_e = 0.01;

t = 0:dt:100;
dt = 0.01; 

i_m = zeros(1,length(t));
n = zeros(1,length(t));
m = zeros(1,length(t));
h = zeros(1,length(t));

V = zeros(1,length(t));
V(1) = (membrane_potential); 

a_n=(0.01*(V(1)+55)) / (1 - exp(-0.1 * (V(1) + 55)));
a_m=(0.1*(V(1)+40)) / (1 - exp(-0.1 * (V(1) + 40)));
a_h=0.07*exp(-0.05*(V(1)+65));

b_n=0.125*exp(-0.0125*(V(1)+65));
b_m=4*exp(-0.0556*(V(1)+65));
b_h=1/(1+exp(-0.1*(V(1)+35)));

n(1)=(a_n*V(1)) / (a_n*V(1) + b_n*V(1));
m(1)=(a_m*V(1)) / (a_m*V(1) + b_m*V(1));
h(1)=(a_h*V(1)) / (a_h*V(1) + b_h*V(1));


%%

for i = 2: length(t)
    
    a_n = (0.01 * (V(i-1)+55)) / (1 - exp(-0.1 * (V(i-1) + 55)));
    a_m = (0.1 * (V(i-1)+40)) / (1 - exp(-0.1 * (V(i-1) + 40)));
    a_h = 0.07 * exp(-0.05*(V(i-1)+65));
    
    b_n = 0.125*exp(-0.0125*(V(i-1)+65));
    b_m = 4*exp(-0.0556*(V(i-1)+65));
    b_h = 1/(1+exp(-0.1*(V(i-1)+35)));
    
    dn = (a_n*(1-n(i-1)) - (b_n*n(i-1)))*dt;
    n(i) = n(i-1)+dn;
    
    dm = (a_m*(1-m(i-1)) - (b_m*m(i-1)))*dt;
    m(i) = m(i-1)+dm;
    
    dh = (a_h*(1-h(i-1)) - (b_h*h(i-1)))*dt;
    h(i) = h(i-1)+dh;
    
    im = (g(1) * (V(i-1) - e(1))) + g(2) * (n(i-1)^4)*(V(i-1)- e(2)) + g(3)*(m(i-1)^3)*h(i-1)*(V(i-1) - e(3));
    dv = (-im + (I_e/A)) * (dt/cm);
    V(i) = V(i-1)+dv;
    
end

%%

figure
subplot(2, 1, 1)
plot(t, V)
ylabel('Membrane Potential (mV)')
xlabel('Time (milliseconds)')
title('Hodgkin−Huxley Model of an Action Potential with 0.01 Microamp Injection Current')

subplot(2, 1, 2)
plot(t, n)
hold on 
plot(t, m)
plot(t, h)
title('Activation Probabilities v. Time')
xlabel('Time (milliseconds)')
ylabel('Probability')
legend('n(t)','m(t)','h(t)')

%%
% <latex>
% 	\item Provide an interesting figure along with an explanation illustrating the extension. (4 pts)
% </latex>

%%

% Top plot: Hodgkin-Huxley Model of an Action Potential. When the membrane potential
% rises above the threshold a spike occurs and the potential decreases until another action potential can be generated. 

% Bottom plot: the activation probabilities of the n, m, and h gating
% variables. Each probability increases and decreases with time. These
% correspond to the membrane potential and activation potential. 

%%
% <latex>
%   \end{enumerate}
% \end{enumerate}
% </latex>


