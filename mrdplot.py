import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import mapvbvd
sys.path.insert(0, '../../mrd-fork/python')
import mrd
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

current_filename = ""
current_header = []
current_streamables_list = []
current_pulse_list = []
current_gradient_list = []
current_acquisition_list = []
current_image_list = []

refresh_plot = 0

fname = 'C:\\Users\\steph\\Desktop\\data\\shurik_all_ischemia_data_081425\\ischemia_27_1\\epsi\\raw.mrd2'
with mrd.BinaryMrdReader(fname) as w:
    current_header = w.read_header()
    current_streamables_list = list(w.read_data())
    current_pulse_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Pulse]
    current_pulse_list.sort(key = lambda x: x.head.pulse_time_stamp_ns)
    current_gradient_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Gradient]
    current_gradient_list.sort(key = lambda x: x.head.gradient_time_stamp_ns)
    current_acquisition_list = [x.value for x in current_streamables_list if type(x.value) == mrd.Acquisition]
    current_acquisition_list.sort(key = lambda x: x.head.acquisition_time_stamp_ns)
plt.figure()
if(current_pulse_list):
    plt.subplot(3,1,1)
    for p in current_pulse_list:
        for q in range(p.amplitude.shape[0]):
            plt.plot(((np.array(range(p.amplitude.shape[1] + 2)) - 1) * p.head.sample_time_ns + p.head.pulse_time_stamp_ns) * \
                     1.0E-9, np.concatenate((np.array([0]), p.amplitude[q, :], np.array([0]))))
    plt.xlabel('time (s)')
    plt.ylabel('pulse amplitude (V)')
if(current_gradient_list):
    plt.subplot(3,1,2)
    for g in current_gradient_list:
        plt.plot((np.array(range(len(g.rl))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
                 1.0E-9, g.rl)
        plt.plot((np.array(range(len(g.ap))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
                 1.0E-9, g.ap)
        plt.plot((np.array(range(len(g.fh))) * g.head.gradient_sample_time_ns + g.head.gradient_time_stamp_ns) * \
                 1.0E-9, g.fh)
    plt.xlabel('time (s)')
    plt.ylabel('pulse amplitude (V)')
if(current_acquisition_list):
    plt.subplot(3,1,3)
    for a in current_acquisition_list:
        plt.plot(np.array(range(len(a.data))) * a.head.sample_time_ns + a.head.acquisition_time_stamp_ns, np.real(a.data[:,0]))
        plt.plot(np.array(range(len(a.data))) * a.head.sample_time_ns + a.head.acquisition_time_stamp_ns, np.imag(a.data[:,0]))
    plt.show()
