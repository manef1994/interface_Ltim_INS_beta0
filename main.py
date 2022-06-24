import antropy as ant
import bisect
from tkinter import *
from tkinter import ttk, messagebox, scrolledtext, filedialog
import tkinter.filedialog
from tkinter.ttk import Progressbar
import numpy as np
import time
import matplotlib
import progress
from scipy.signal import butter, iirnotch, lfilter, filtfilt
import matplotlib.pyplot as plt
import neurokit2 as nk
from shapely.geometry import LineString
from mne.io import read_raw_edf

matplotlib.pyplot.figure
matplotlib.lines
matplotlib.lines.Line2D

def butter_lowpass(cutoff, fs, order=2):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_highpass(cutoff, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def remove_baseline_wander(data, sample_rate, cutoff=0.05):
    return filter_signal(data=data, cutoff=cutoff, sample_rate=sample_rate, filtertype='notch')

def butter_bandpass(lowcut, highcut, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def filter_signal(data, cutoff, sample_rate, order=2, filtertype='lowpass', return_top=False):
    if filtertype.lower() == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
    cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    elif filtertype.lower() == 'notch':
        b, a = iirnotch(cutoff, Q=0.005, fs=sample_rate)
    else:
        raise ValueError('filtertype: %s is unknown, available are: \
    lowpass, highpass, bandpass, and notch' % filtertype)

    filtered_data = filtfilt(b, a, data)

    if return_top:
        return np.clip(filtered_data, a_min=0, a_max=None)
    else:
        return filtered_data

def std_compute(xx):
    i = 0
    j = 6
    STD = []

    while (True):
        go_in = []

        go_in = xx[i:j]
        std = np.std(go_in)
        STD.append(std)

        j = j + 6
        i = i + 6

        if j > len(xx):
            break
    return STD

def Browse():
    plt.close(1)
    plt.close(2)
    plt.close(3)
    plt.close(4)
    plt.close(5)
    plt.close(6)
    plt.close(7)
    plt.close(8)
    plt.close(9)
    plt.close(10)

    sText.delete(1.0, END)
    progress['value'] = 0

    global path_signal
    global EDF

    # ======= open file selector
    main_win.sourceFile = tkinter.filedialog.askopenfilename(parent=main_win, initialdir="E:\\data\\siena-scalp-eeg-database-1.0.0\\",
                                                             title='Please select a directory',)
    path = main_win.sourceFile
    path = path.replace('/', '\\')
    # ======= print path

    if path.endswith('.csv') or path.endswith('.edf'):
        if path.endswith('.edf'):

            e1.config(state='disabled')
        if path.endswith('.csv'):
            e1.config(state='normal')

        path_browse.config(text=path)
        path_signal = str(path)

        sText.insert(END, 'signal path equale to: ' + str(path) + '\n')
        return path_browse.config(text=path_signal)

    else:
        messagebox.showinfo("warninig file missing",
                            "you must select a CSV or EDF file, please recheck the file selected")

def Read_signal(path):
    global fs
    signal_input = []

    if path.endswith('.csv'):
        sText.insert(END, 'reading the signal... \n')
        time.sleep(1)
        signal_input = np.loadtxt(path)
        fs = int(e1.get())

    elif path.endswith('.edf'):

        edf = read_raw_edf(path, preload=False, stim_channel=None, verbose=False)
        xx = edf.ch_names
        index = xx.index("2")
        fs = edf.info['sfreq']
        fs = int(fs)
        signal_input = edf[index]
        signal = signal_input[0]
        print("FS out of the loop equals to: ", fs)

    signal_input = []
    signal_input = signal[0]

    sText.insert(END, 'ECG signal is in the channel: ' + str((index) + 1) + '\n')
    sText.insert(END, 'lenght of the ECG signal equal to ' + (str(((len(signal_input)) / fs) / 60)) + ' minutes ' + '\n')
    time.sleep(1)

    print('signal length is:', (len(signal_input)/fs)/60)

    # == reading the exact time of the seizure tapped by the user

    seizure_time = str(seizure_start.get())

    res = [float(value) for value in seizure_time.split(':')]

    print("res \t", res)
    print("res \t", int(res[0]))
    print("res \t", int(res[1]))
    print("res \t", int(res[2]))

    seizure = (fs * 60 * 60 * int(res[0])) + (fs * 60 * int(res[1])) + (fs * int(res[2]))
    #seizure = (fs * 60 * 60 * int(res[0])) + (fs * 60 * int(res[1])) + (fs * int(res[2]))

    pre_ictal = seizure - (fs * 60 * 60 * 1)
    end = seizure + (fs * 60 * 10)

    print('fs:\t', fs)
    print('pre-ictal:\t', pre_ictal)
    print('seizure:\t', seizure)
    print('end:\t', end)

    signal_input = signal_input[pre_ictal:end]

    print("signal length\t", len(signal_input))

    return signal_input

def r_algorithm(signal, fs):

    detection_algo = combo_r_peak.get()
    print('detectrion algorithm selected\t', detection_algo)
    r_peaks = []

    if (detection_algo == 'pantompkins1985'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="pantompkins1985")

    if (detection_algo == 'hamilton2002'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="hamilton2002")

    if (detection_algo == 'christov2004'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="christov2004")

    if (detection_algo == 'elgendi2010'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="elgendi2010")

    if (detection_algo == 'christov2004'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="christov2004")

    if (detection_algo == 'rodrigues2021'):
        r_peaks = nk.ecg_peaks(signal, sampling_rate=fs, method="rodrigues2021")
        print("peaks length before return\t", len(r_peaks))

    xx = r_peaks[1]
    r_peaks = xx["ECG_R_Peaks"]
    print("first peak\t", r_peaks[0])
    print("final peak\t", r_peaks[-1])

    # pl = [0] * len(r_peaks)
    #
    # fig, axs = plt.subplots()
    # axs.plot(signal_input)
    # axs.plot(r_peaks, pl, 'ro')
    # axs.set_xlabel('Time')
    # axs.set_ylabel('Amplitude')
    # axs.grid(True)
    # plt.show()
    #
    # print("peaks length before return2\t", len(r_peaks))
    return r_peaks

def Action(path):
    global fs, seizures_start_f, final_peaks, signal_input
    #plt.close('all')
    first_filter = comboExample.get()
    seconde_filter = filter_two.get()
    third_filter = filter_three.get()

    detection_algo = combo_r_peak.get()

    cuttoff_int = float(cutoff.get())
    cuttoff_int_two = float(cuttoff_two.get())
    cuttoff_int_three = float(cuttoff_three.get())

    thre_ApEn = thresh_ApEn.get()
    thre_NN50 = thresh_NN50.get()


    if (cuttoff_int_two != 0) and (cuttoff_int != 0) and (cuttoff_int_three != 0):

        sText.insert(END, 'Reading ECG file.....\n')
        print('reading ECG signal')
        signal_input = Read_signal(path)
        print('fs equals to after signal reading: ', fs)
        print('signal intput length', len(signal_input))
        time.sleep(1)

        sText.insert(END, 'first filtering: ' + str(first_filter) + " " + " with a cutoff equal to: " + str(cuttoff_int) + '\n')
        signal_output = filter_signal(data=signal_input, cutoff=cuttoff_int, sample_rate=fs, order=2, filtertype=first_filter)

        progress['value'] = 10
        window.update_idletasks()
        time.sleep(1)

        sText.insert(END, 'second filtering: ' + str(seconde_filter) + " " + " with a cutoff equal to: " + str(cuttoff_int_two) + '\n')
        signal_output = filter_signal(data=signal_output, cutoff=cuttoff_int_two, sample_rate=fs, order=2,filtertype=seconde_filter)

        progress['value'] = 20
        window.update_idletasks()
        time.sleep(1)

        sText.insert(END, 'third filtering: ' + str(third_filter) + " " + " with a cutoff equal to: " + str(cuttoff_int_three) + '\n')
        signal_output = filter_signal(data=signal_output, cutoff=cuttoff_int_three, sample_rate=fs, order=2, filtertype=third_filter)

        progress['value'] = 25
        window.update_idletasks()
        time.sleep(1)

        #######################################################################################################
        # == features extraction
        sText.insert(END, 'features extraction....\n')
        print('features extraction....\n')
        time.sleep(1)

        ApEN, NN50 = features_EX(signal_output, fs)
        print('featyures extraction done')

        progress['value'] = 50
        window.update_idletasks()
        time.sleep(1)

        # == computing the STD curves

        STD_app = std_compute(ApEN)
        STD_NN50 = std_compute(NN50)

        progress['value'] = 90
        window.update_idletasks()
        time.sleep(1)

        thresh_AP = thre_ApEn
        thresh_nn = thre_NN50

        print('ApEn threshold equals to:\t', thresh_AP)
        print('NN50 threshold equals to:\t', thresh_nn)
        # =============================================================================================================
        # == finding the intersections of the threshold value computed

        first = 0
        tt = [first]
        i = 1
        while i <= len(ApEN) - 1:
            first = first + 10
            tt.append(first)
            i = i + 1

        first = 0
        tt_healthy = [first]
        i = 1
        while i <= len(STD_app) - 1:
            first = first + 60
            tt_healthy.append(first)
            i = i + 1

        first = 0
        tt_healthy = [first]
        i = 1
        while i <= len(STD_app) - 1:
            first = first + 1
            tt_healthy.append(first)
            i = i + 1

        # ____ woerking on the approximate thresholding
        x = np.zeros(len(STD_app), dtype=float, order='C')
        xy = np.zeros(len(ApEN), dtype=float, order='C')

        arr = []
        arr1 = []

        for i in range(len(STD_app)):
            arr.append(thresh_AP)
            arr1.append(0)

        first_line = LineString(np.column_stack((tt_healthy, STD_app)))
        second_line = LineString(np.column_stack((tt_healthy, arr)))
        intersection = first_line.intersection(second_line)
        # x, y = LineString(intersection).xy
        # print(" x \t", x)

        arr_n = []
        arr1_n = []

        for i in range(len(STD_NN50)):
            arr_n.append(thresh_nn)
            arr1_n.append(0)

        first_line_n = LineString(np.column_stack((tt_healthy, STD_NN50)))
        second_line_n = LineString(np.column_stack((tt_healthy, arr_n)))
        intersection_n = first_line_n.intersection(second_line_n)
        # x_n, y_n = LineString(intersection_n).xy
        # print((x_n))

        round_app = []
        round_NN50 = []

        # =============================================================================================================
        # == plotting the results

        progress['value'] = 100
        window.update_idletasks()
        time.sleep(1)

        fig, axs = plt.subplots(2, 1)
        axs[0].plot(tt, ApEN, label="approximate entropy of the input signal", marker='o')
        axs[1].plot(tt_healthy, STD_app, label="STD of approximate entropy", marker='o')
        axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
        axs[1].axvline(x=len(STD_app) - 10, color='red', linestyle='--')

        axs[0].set_title('Subject curves', fontsize=24, y=1)

        axs[1].axhline(y=thresh_AP, color='blue', linestyle='--')

        axs[0].set_xlabel('sample per min')
        axs[0].set_ylabel('entropy value')
        axs[1].set_xlabel('sample per min')
        axs[1].set_ylabel('entropy value')
        axs[0].grid(True)
        axs[1].grid(True)
        axs[0].legend()
        axs[1].legend()

        ymin = 0.1
        ymax = 1.3
        axs[0].set_ylim([ymin, ymax])

        ymin_1 = 0
        ymax_1 = 0.25
        axs[1].set_ylim([ymin_1, ymax_1])

        if intersection.geom_type == 'MultiPoint':
            axs[1].plot(*LineString(intersection).xy, 'o')
        elif intersection.geom_type == 'Point':
            axs[1].plot(*intersection.xy, 'o')

        # == plotting NN50
        fig, axs = plt.subplots(2, 1)
        axs[0].plot(tt_healthy, STD_app, label="STD of ApEn feature", marker='o')
        axs[1].plot(tt_healthy, STD_NN50, label="STD of NN50 feature", marker='o')

        axs[0].set_title('Subject curves ', fontsize=24, y=1)

        axs[0].axvline(x=len(STD_app) - 10, color='red', linestyle='--')
        axs[1].axvline(x=len(STD_NN50) - 10, color='red', linestyle='--')

        axs[0].set_xlabel('sample per min')
        axs[0].set_ylabel('entropy value')
        axs[1].set_xlabel('sample per min')
        axs[1].set_ylabel('entropy value')
        axs[0].grid(True)
        axs[1].grid(True)
        axs[0].legend()
        axs[1].legend()

        axs[0].axhline(y=thresh_AP, color='blue', linestyle='--')
        axs[1].axhline(y=thresh_nn, color='blue', linestyle='--')

        ymin_1 = 0;
        ymax_1 = 0.25
        axs[0].set_ylim([ymin_1, ymax_1])

        ymin_1 = 0;
        ymax_1 = 14
        axs[1].set_ylim([ymin_1, ymax_1])

        if intersection.geom_type == 'MultiPoint':
            axs[0].plot(*LineString(intersection).xy, 'o')

        elif intersection.geom_type == 'Point':
            axs[0].plot(*intersection.xy, 'o')

        if intersection_n.geom_type == 'MultiPoint':
            axs[1].plot(*LineString(intersection_n).xy, 'o')

        elif intersection_n.geom_type == 'Point':
            axs[1].plot(*intersection_n.xy, 'o')

        # == plotting NN50
        fig, axs = plt.subplots(2, 1)
        axs[0].plot(tt, NN50, label=" NN50 of the input signal", marker='o')
        axs[1].plot(tt_healthy, STD_NN50, label="STD of NN50 feature", marker='o')
        axs[0].set_title('Subject curves ', fontsize=24, y=1)

        axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
        axs[1].axvline(x=len(STD_NN50) - 10, color='red', linestyle='--')

        axs[0].set_xlabel('sample en seconde')
        axs[0].set_ylabel('valeur de entropie')
        axs[1].set_xlabel('sample en seconde')
        axs[1].set_ylabel('valeur de entropie')
        axs[0].grid(True)
        axs[1].grid(True)
        axs[0].legend()
        axs[1].legend()

        axs[1].axhline(y=thresh_nn, color='blue', linestyle='--')

        ymin = 70;
        ymax = 275
        axs[0].set_ylim([ymin, ymax])

        ymin_1 = 0;
        ymax_1 = 14
        axs[1].set_ylim([ymin_1, ymax_1])

        if intersection_n.geom_type == 'MultiPoint':
            axs[1].plot(*LineString(intersection_n).xy, 'o')
        elif intersection_n.geom_type == 'Point':
            axs[1].plot(*intersection_n.xy, 'o')

        sText.insert(END, "--- %s seconds ---" % (time.time() - start_time))

        plt.show()

        # ==  formatters initialisation
        formatter = matplotlib.ticker.FuncFormatter(lambda ms, x: time.strftime('%M:%S', time.gmtime(ms / fs)))

    else:
        messagebox.showinfo("warninig value missing", " sample rate and cutoff value must be greater than 0 !")

def features_EX(singal_ouput, fs):


    print('how many minute: ', int(len(singal_ouput) / ((fs * 6) * 10)))

    # == R peaks detection
    peaks = r_algorithm(singal_ouput, fs)

    NN50, ApEn, = ([] for i in range(2))
    start = 0
    end = 120 * fs


    print('peaks length:\t', len(peaks))

    while True:

        go = bisect.bisect_left(peaks, start)
        out = bisect.bisect_left(peaks, end)

        RRi = []
        peaks_in = []
        ff = 0
        peaks_in = peaks[go:out]

        for i in range(len(peaks_in) - 1):
            new = peaks_in[i + 1] - peaks_in[i]
            new = (new / fs)
            RRi.append(new)
            if (new > 0.05):
                ff = ff + 1

        NN50.append(ff)
        ApEn.append(ant.app_entropy(RRi))

        start = start + (10 * fs)
        end = end + (10 * fs)
        i += 1

        if (end > peaks[-1]):
            break

    return ApEn, NN50

    formatter = matplotlib.ticker.FuncFormatter(lambda ms, x: time.strftime('%H:%M:%S', time.gmtime((ms / fs)/6)))

signal_input = []
start_time = time.time()
heart_rate_bmp = []
HRV, mean_hr, max_hr, min_hr, average_diff = ([] for i in range(5))
features_num = 0
feature_nn50 = 0
sdd = 0
fs = 0
final_peak = []
ApEn = []
NN50 = []

row = 0

# create new window object NP: all the interface code must be between Tk() and mainloop
EDF = False
window = Tk()
# define the window standard size
window.geometry('1300x600')
window.title("ECG manipulation ")
# define labels

# initiate tinker and hide window
main_win = tkinter.Tk()
main_win.withdraw()

main_win.overrideredirect(True)
main_win.geometry('0x0+0+0')

main_win.deiconify()
main_win.lift()
main_win.focus_force()

# ======= take the path of the file to use
l1 = Label(window, text="EDF/CSV file to select")
l1.grid(row=row, column=0)


button = Button(text="Browse", fg="red", command=Browse, height=2, width=30)
button.grid(row=row, column=1)

path_browse = Label(window, text="the path of your file is")
path_browse.grid(row=row, column=2)
row += 1

space = Label(window, text="      ")
space.grid(row=row, column=0)
row += 1
# ======= take the frequency of the input signal
l2 = Label(window, text="frequency sample in Hz")
l2.grid(row=row, column=0)

e1 = Entry(window, width="15")
e1.grid(row=row, column=1)
e1.insert(END, '0')
row += 1

space = Label(window, text="      ")
space.grid(row=row, column=0)
row += 1

# ======= select the first filter to be used
fil = Label(window, text="the filter to use")
fil.grid(row=row, column=0)

# Combobox of the filter selection
fontExample = ("Courier", 16, "bold")
comboExample = ttk.Combobox(window, values=["highpass", "lowpass", "notch", ])
comboExample.grid(row=row, column=1)
comboExample.current(2)

# cutoff
fil1 = Label(window, text="the frequency to use by the filter")
fil1.grid(row=row, column=2)
cutoff = Entry(window, width="35")
cutoff.grid(row=row, column=3)
cutoff.insert(END, '60')
row += 1

# ======= select the second filter to be used
fil = Label(window, text="the filter to use")
fil.grid(row=row, column=0)

# Combobox of the filter selection
filter_two = ttk.Combobox(window, values=["highpass", "lowpass", "notch", ])
filter_two.grid(row=row, column=1)
filter_two.current(0)

# cutoff
space = Label(window, text="      ")
space.grid(row=row, column=0)
fil1 = Label(window, text="the frequency to use by the filter")
fil1.grid(row=row, column=2)
cuttoff_two = Entry(window, width="35")
cuttoff_two.grid(row=row, column=3)
cuttoff_two.insert(END, '10')
row += 1
# ======= select the third filter to be used
fil = Label(window, text="the Filter to use")
fil.grid(row=row, column=0)

# Combobox of the filter selection
filter_three = ttk.Combobox(window, values=["highpass", "lowpass", "notch", ])
filter_three.grid(row=row, column=1)
filter_three.current(1)

# cutoff
space = Label(window, text="      ")
space.grid(row=row, column=0)
fil1 = Label(window, text="the frequency to use by the filter")
fil1.grid(row=row, column=2)
cuttoff_three = Entry(window, width="35")
cuttoff_three.grid(row=row, column=3)
cuttoff_three.insert(END, '20')

row += 1
space = Label(window, text="      ")
space.grid(row=row, column=0)
row += 1

# ======= select the r peak algorithm to be used
R_pe = Label(window, text="R peak algorithm to use")
R_pe.grid(row=row, column=0)

# Combobox
combo_r_peak = ttk.Combobox(window, values=["pantompkins1985", "hamilton2002", "christov2004", "elgendi2010", "engzee2012", "rodrigues2021"])

combo_r_peak.grid(row=row, column=1)
combo_r_peak.current(5)
row += 1

#== drawing the seizures line

space = Label(window, text="      ")
space1 = Label(window, text="      ")
space.grid(row=row, column=0)
row+=1
space1.grid(row=row, column=0)
sei = Label(window, text="The seizure start at ")
sei.grid(row=row, column=0)
seizure_start = Entry(window, width="35")
seizure_start.grid(row=row, column=1)
seizure_start.insert(END, 'h:m:s')
row += 1

# ======= reading threshold values
# = ApEn threshold value
space = Label(window, text="      ")
space.grid(row=row, column=0)
row += 1

space = Label(window, text="      ")
space1 = Label(window, text="      ")
space.grid(row=row, column=0)
row+=1
space1.grid(row=row, column=0)
sei = Label(window, text="The threshold value for ApEn feature")
sei.grid(row=row, column=0)
thresh_ApEn = Entry(window, width="35")
thresh_ApEn.grid(row=row, column=1)
thresh_ApEn.insert(END, '0')
# row += 1

# # = NN50 threshold value
# space = Label(window, text="      ")
# space.grid(row=row, column=0)
# row += 3
#
space = Label(window, text="      ")
space1 = Label(window, text="      ")
space.grid(row=row, column=0)
row+=1
space1.grid(row=row, column=0)
sei = Label(window, text="The threshold value for NN50 feature")
sei.grid(row=row, column=0)
thresh_NN50 = Entry(window, width="35")
thresh_NN50.grid(row=row, column=1)
thresh_NN50.insert(END, '0')
row += 1



# ======= start the work
# space = Label(window, text="      ")
# space.grid(row=row, column=0)
# row += 1

button_action = Button(text="Apply", fg="red", height=3, width=60, command=lambda: Action(path=path_signal))
button_action.grid(row=row, column=1)

sText = scrolledtext.ScrolledText(window, height=10)
sText.insert(INSERT, "steps Progess...\n")
sText.insert(END, "")
sText.grid(column=2, row=row, sticky='WE', columnspan=3)
row += 1

space = Label(window, text="      ")
space.grid(row=row, column=0)
row += 1

Pro = Label(window, text="Software Progress... ")
Pro.grid(row=row, column=0)

progress = Progressbar(window, orient=HORIZONTAL, length=1075, mode='determinate')
progress.grid(row=row, column=1, columnspan=4)
row += 1

space = Label(window, text="      ")
space.grid(row=row, column=0)
copy_right = Label(window,
                   text=" \N{COPYRIGHT SIGN} \N{TRADE MARK SIGN} \N{REGISTERED SIGN} all right reserved: LTIM by Manef Ben Mbarek")
copy_right.grid(row=row, column=3, columnspan=3)
copy_right.config(font=("Courier", 10))

window.mainloop()

__author__ = "Manef Ben mbarek"
__copyright__ = "Copyright (C) 2020 Manef"
__license__ = "Public Domain"
__version__ = "1.0"
