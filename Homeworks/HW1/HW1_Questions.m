    %%
% <latex>
% \title{BE 521: Homework 1 \\{\normalsize Exploring Neural Signals} \\{\normalsize Spring 2021}}
% \author{33 points}
% \date{Due: Tuesday 2/2/2021 10 PM}
% \maketitle
% \textbf{Objective:} Working with the IEEG Portal to explore different Neural signals
% </latex>

%% 
% <latex>
% \begin{center}
% \author{Saif Khawaja \\
%   \normalsize Collaborators: Raveen K\\}
% \end{center}
% </latex>

%%
% <latex>
% \section{Seizure Activity (16 pts)} 
% The dataset \texttt{I521\_A0001\_D002} contains an example of human intracranial EEG (iEEG) data displaying seizure activity. It is recorded from a single channel (2 electrode contacts) implanted in the hippocampus of a patient with temporal lobe epilepsy being evaluated for surgery. In these patients, brain tissue where seizures are seen is often resected. You will do multiple comparisons with this iEEG data and the unit activity that you worked with in Homework 0 \texttt{(I521\_A0001\_D001)}. You will have to refer to that homework and/or dataset for these questions.
% \begin{enumerate}
%  \item Retrieve the dataset in MATLAB using the IEEGToolbox and generate a \emph{session} variable as before (No need to report the output this time). What is the sampling rate of this data? What is the maximum frequency of the signal content that we can resolve? (2 pts)
% </latex>

%% 

% The sampling rate is 200 Hz. We can resolve a maximum frequency of 200 Hz
% / 2 = 100 Hz.

session = IEEGSession('I521_A0001_D002', 'saifkhawaja', 'sai_ieeglogin.bin')

session.data

Hz = session.data.sampleRate

%%
% <latex>
%  \item How does the duration of this recording compare with the recording from HW0 \texttt{(I521\_A0001\_D001)}? (2 pts)
% </latex>

% Obtain the duration of the recording for channel 1 (retuyrned in microseconds).

durationInUSec = session.data(1).rawChannels(1).get_tsdetails.getDuration;

durationInSec = durationInUSec / 1e6

% durationInSec = 644.9950 seconds, which contrasts with 10s recorded in HW0.


%%
% <latex>
%  \item Using the time-series visualization functionality of the IEEG Portal, provide a screenshot of the first 500 ms of data from this recording. (2 pts)
% </latex>

% \includegraphics[scale=0.3]{/Users/saif/Documents/GitHub/braincomputerinterfaces/Homeworks/HW1/signal.png}\\

%%
% <latex>
%  \item Compare the activity in this sample with the data from HW0.  What differences do you notice in the amplitude and frequency characteristics? (2 pts)
% </latex>

% Obtain HW1 data initially: 

nr = ceil((session.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session.data.sampleRate); 
allData = session.data.getvalues(1:nr,1);
sr2 = 200;

ms1 = allData(1:0.5*sr2); 
peak1 = max(ms1)

% Initialize session and pull HW0 data: 

session0 = IEEGSession('I521_A0001_D001', 'saifkhawaja', 'sai_ieeglogin.bin')

session0.data

sr1 = 32051;
nr = ceil((session0.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session.data.sr1); 
allData = session0.data.getvalues(1:nr,1);
sr1 = 200;

ms0 = allData(1:0.5*sr1); 
peak0 = max(ms1)

A_D001 = mean(ms0) 
A_D002 = mean(ms1)

%%
% <latex>
%  \item The unit activity sample in \texttt{(I521\_A0001\_D001)} was high-pass filtered to remove low-frequency content. Assume that the seizure activity in \texttt{(I521\_A0001\_D002)} has not been high-pass filtered. Given that the power of a frequency band scales roughly as $1/f$, how might these differences in preprocessing contribute to the differences you noted in the previous question? (There is no need to get into specific calculations here. We just want general ideas.) (3 pts)
% </latex>

% <latex>
%  \item The removal of noise through the high-pass filtering in pre-processing cleaning of the data results in more prominent striations in the data. There are more oscillations in D001 than in D002, which would likely be resolved if the same pre-processing was implemented in D002.
% Similarly, if we had preprovessed D002, we would observe flatter
% oscillations and lower amplitudes because, as reference in the 1/f in
% question, the power scales exponentially as frequency shrinks as the
% amplitude increases. This means that the waves on the reading have lower
% power and would flatten out if removed. Compared to D001, this would be flatter and if we did not remove the lower frequencies in D001 we would observe larger amplitude and sharper disparations in spikes as the signal would have a larger amount of content in lower frequencies.  
% </latex>

%%
% <latex>
%  \item Two common methods of human iEEG are known as electrocorticography (ECoG) and stereoelectroencephalography (SEEG). For either of these paradigms (please indicate which you choose), find and report at least two of the following electrode characteristics: shape, material, size. Please note that exact numbers aren't required, and please cite any sources used. (3 pts)
% </latex>

% Electrocorticography (ECoG) places electrodes directly on exposed
% surfaces of the brain to record cerebral cortex activity. ECoG electrode
% arrays are typically organized as sixteen strile stainless steel carbon
% tipped electrodes made of platinum, platinum-indium alloys or gold. Grid
% electrodesa are widely used and have between 4 and 256 contacts.
% Spacing is normally 1cm and individual electrodes are ~5mm in diameter.
% References:
% Mesgarani, N; Chang, EF (2012). "Selective cortical representation of attended speaker in multi-talker speech perception". Nature. 485 (7397): 233–6. Bibcode:2012Natur.485..233M. doi:10.1038/nature11020.
% Schuh, L; Drury, I (1996). "Intraoperative electrocorticography and direct cortical electrical stimulation". Seminars in Anesthesia. 16: 46–55. doi:10.1016/s0277-0326(97)80007-4

%%
% <latex>
%  \item What is a local field potential? How might the  characteristics of human iEEG electrodes cause them to record local field potentials as opposed to multiunit activity, which was the signal featured in HW0 as recorded from 40 micron Pt-Ir microwire electrodes? (2 pts)
% </latex>

% A local field potential is an extracellular signal generated from ion
% concentration imbalances outside of cells. The potentials are generated
% by the voltage from charge separation across this ion concentration
% difference (relative excitatory and inhibitory dendric potentials).
%
% iEEG electrodes can record local field potentials because of their
% comparitively large size (40 micron vs. 5mm diameter). This means the
% signals recorded and read are a summation of multiple neuron potentials. The 40 micron Pt-Ir electrodes pick up signals from multiunit activity because of their size and the much higher sampling rate required for recording the minutia in adjustments.

% Reference: Stacey WC, Kellis S, Greger B, et al. Potential for unreliable
% interpretation of EEG recorded with microelectrodes. Epilepsia.
% 2013;54(8):1391-1401. doi:10.1111/epi.12202 

%%
% <latex>
% \end{enumerate}
% </latex>

%%
% <latex>
% \section{Evoked Potentials (17 pts)} 
% The data in \texttt{I521\_A0001\_D003} contains an example of a very common type of experiment and neuronal signal, the evoked potential (EP). The data show the response of the whisker barrel cortex region of rat brain to an air puff stimulation of the whiskers. The \texttt{stim} channel shows the stimulation pattern, where the falling edge of the stimulus indicates the start of the air puff, and the rising edge indicates the end. The \texttt{ep} channel shows the corresponding evoked potential. 
% Once again, play around with the data on the IEEG Portal, in particular paying attention to the effects of stimulation on EPs. You should observe the data with window widths of 60 secs as well as 1 sec. Again, be sure to explore the signal gain to get a more accurate picture. Finally, get a sense for how long the trials are (a constant duration) and how long the entire set of stimuli and responses are.
% </latex>

%%
% <latex>
% \begin{enumerate}
%  \item Based on your observations, should we use all of the data or omit some of it? (There's no right answer, here, just make your case either way in a few sentences.) (2 pts)
% </latex>

% While there is a large similarity in most of the readings that there is a
% noticeable peak, many of the recorded signals also have a large trough
% (and sometimes multiple trouphs). Moreover, some trials had multiple peaks.
% These indicate there could be additional stimuli (lights, noises, and
% other sensory inputs) that triggered a response in the rat, and therefore 
% these recordings should be removed as they do not isolate the signal we are 
% looking for. 
%%
% <latex>
%  \item Retrieve the \texttt{ep} and \texttt{stim} channel data in MATLAB. What is the average latency (in ms) of the peak response to the stimulus onset over all trials? (Assume stimuli occurs at exactly 1 second intervals)(3 pts)
% </latex>

% Initialize session for dataset 003 and pull information

session3 = IEEGSession('I521_A0001_D003', 'saifkhawaja', 'sai_ieeglogin.bin');

session3.data(1)

nr = ceil((session3.data.rawChannels(1).get_tsdetails.getEndTime)/1e6*session3.data.sampleRate);
data3 = session3.data.getvalues(1:nr, 1:2);

ep = data3(:,1);
stim = data3(:,2);

% Obtain durtation and sample rate

durationInUSec = session3.data(1).rawChannels(1).get_tsdetails.getDuration; 
durationInSec = durationInUSec / 1e6;
sr = session3.data.sampleRate;

%%
% <latex>
%  \item In neuroscience, we often need to isolate a small neural signal buried under an appreciable amount of noise.  One technique to accomplish this is called the spike triggered average, sometimes called signal averaging. This technique assumes that the neural response to a repetitive stimulus is constant (or nearly so), while the noise fluctuates from trial to trial - therefore averaging the evoked response over many trials will isolate the signal and average out the noise.
%  Construct a spike triggered average plot for the data in \texttt{I521\_A0001\_D003}.  Plot the average EP in red.  Using the commands \texttt{hold on} and \texttt{hold off} as well as \texttt{errorbar} and \texttt{plot}, overlay error bars at each time point on the plot to indicate the standard deviation of the responses at any given time point.  Plot the standard deviation error bars in gray (RGB value: [0.7 0.7 0.7]). Make sure to give a proper legend along with your labels. (4 pts)
% </latex>

%%
% <latex>
%  \item 
%   \begin{enumerate}
% 	\item We often want to get a sense for the amplitude of the noise in a single trial. Propose a method to do this (there are a few reasonably simple methods, so no need to get too complicated). Note: do not assume that the signal averaged EP is the ``true'' signal and just subtract it from that of each trial, because whatever method you propose should be able to work on the signal from a single trial or from the average of the trials. (4 pts)
% </latex>

%%
% <latex>
% 	\item Show with a few of the EPs (plots and/or otherwise) that your method gives reasonable results. (1 pt)
% </latex>

%%
% <latex>
% 	\item 
%     \begin{enumerate}
%         \item Apply your method on each individual trial and report the mean noise amplitude across all trials. (1 pt)
% </latex>

%%
% <latex>
%         \item Apply your method on the signal averaged EP and report its noise. (1 pt)
% </latex>

%%
% <latex>
% 	    \item Do these two values make sense? Explain. (1 pt)
% </latex>

%%
% <latex>
%     \end{enumerate}
%   \end{enumerate}
% \end{enumerate}
% </latex>

