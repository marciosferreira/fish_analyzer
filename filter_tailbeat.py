import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description='filtering results from fish analyzer software')

# General analysis arguments
parser.add_argument('-i', '--indir',     action="store",         dest='indir',
                    help='Input directory with csv file',                                 default=False,      required=True)
parser.add_argument('-o', '--outdir',    action="store",         dest='outdir',
                    help='Output directory to save reports',                              default=False,      required=True)
parser.add_argument('-s', '--start',     action="store",         dest='start_frame',
                    help='the frame to start analyses',                      default=0,   required=False)
parser.add_argument('-l', '--last',     action="store",         dest='last_frame',
                    help='the last frame to analyze',                       default=-1,       required=False)
parser.add_argument('-f', '--fps',       action="store",         dest='fps',
                    help='Frames per second of the data',                               default=0.0,        required=True,   type=int)


args = parser.parse_args()

in_dir = args.indir  # "C:/Users/marci/Downloads/dugale.csv"
out_dir = args.outdir  # = "C:/Users/marci/Downloads/"
fps = args.fps
start_frame = int(args.start_frame)
finish_frame = int(args.last_frame)

print("")
print("")
print("####################################")
print("Fish Analyzer correction script")
print("in dir = " + in_dir)
print("output dir = " + out_dir)
print("fps = " + str(fps))
print("start frame = " + str(start_frame) +
      ' (0 is the first frame, and is the default value)')
print("finalframe = " + str(finish_frame) +
      ' (-1 means lat frame, and is the default value)')
print("Look at the out_dir for a pdf and a png file with results and reports")
print("####################################")
print("")
print("")


def peaks(angle_list, fps, out_dir):

    mean_points = 0.3

    dataset = pd.DataFrame(angle_list)

    dataset.columns = ['angle']

    dataset = dataset['angle']

    totaltime = dataset.count() / fps
    total_frames = int(dataset.count())
    import matplotlib
    # matplotlib.use('tkagg')
    import matplotlib.pyplot as plt
    import statistics

    # plt.title("Caudal beating rate") #The title of our plot
    # plt.plot(dataset.angle) #Draw the plot object

    import numpy as np
    import math

    # Calculate moving average with 0.75s in both directions, then append do dataset
    hrw = mean_points  # One-sided window size, as proportion of the sampling frequency

    mov_avg = dataset.rolling(
        int(hrw*fps)).mean()  # Calculate moving average

    # Impute where moving average function returns NaN, which is the beginning of the signal where x hrw

    avg = (np.mean(abs(dataset)))
    avg_hr = 0
    mov_avg = [avg_hr if 1 == 1 else x for x in mov_avg]

    window1 = []
    peaklist1 = []
    fps_skip = []
    window2 = []
    peaklist2 = []
    listpos = 0  # We use a counter to move over the different data columns
    comum = 0
    for datapoint in dataset:

        if (datapoint > 0):

            window1.append(datapoint)

        else:
            if len(window1) >= 2:
                #maximum = max(window1)
                # Notate the position on the X-axis
                beatposition = listpos - \
                    len(window1) + (window1.index(max(window1)))
                if max(window1) > (avg*0.3):
                    if comum > max(window1):
                        pass
                    elif comum > 0 and comum < max(window1):
                        peaklist1[-1] = beatposition
                        comum = max(window1)
                    else:
                        peaklist1.append(beatposition)
                        comum = max(window1)
                        fps_skip.append(len(window1))

                window1 = []

            else:
                window1 = []  # Clear marked ROI

        if (datapoint < 0):
            window2.append(datapoint)

        else:
            if len(window2) >= 2:
                #minimum = min(window2)
                # Notate the position on the X-axis
                beatposition = listpos - \
                    len(window2) + (window2.index(min(window2)))
                if min(window2) < ((avg*-1)*0.3):
                    if comum < min(window2):
                        pass
                    elif comum < 0 and comum > min(window2):
                        peaklist2[-1] = beatposition
                        comum = max(window2)
                    else:
                        peaklist2.append(beatposition)
                        comum = min(window2)
                        fps_skip.append(len(window2))

                window2 = []

            else:
                window2 = []  # Clear marked ROI

        listpos += 1

    fps_skip_total = np.mean(fps_skip)

    if fps_skip_total > 4:

        skip_state = ' The frame rate applyed seems ok, with about ' + \
            str(round(fps_skip_total, 2)) + ' in each peak'

    else:

        skip_state = ' Be carefull: the frame rate applyed seems low, with about ' + \
            str(round(fps_skip_total, 2)) + ' in each peak'

    def cvariation(x): return np.std(x, ddof=1) / np.mean(x) * 100

    # Get the y-value of all peaks for plotting purposes
    ybeat1 = [dataset[x] for x in peaklist1]
    cv_ybeat1 = cvariation(ybeat1)

    # Get the y-value of all peaks for plotting purposes
    ybeat2 = [dataset[x] for x in peaklist2]

    cv_ybeat2 = cvariation(ybeat2)
    cv_ybeat3 = (cv_ybeat1 + abs(cv_ybeat2))/2

    peaklist_t = peaklist1 + peaklist2
    ybeat_t = ybeat1 + ybeat2

    amp = np.mean(ybeat1)
    soma = np.sum(ybeat1)
    ampx = round(amp, 2) * 2
    ms_dist1 = 0
    RR_list = []
    cnt = 0
    while (cnt < (len(peaklist1)-1)):
        # Calculate distance between beats in # of samples
        RR_interval = (peaklist1[cnt+1] - peaklist1[cnt])
        # Convert sample distances to ms distances
        ms_dist = ((RR_interval / fps) * 1000.0)
        RR_list.append(ms_dist)  # Append to list
        cnt += 1

    amp2 = np.mean(ybeat2)
    soma2 = np.sum(ybeat2)
    ampy = round(amp2, 2) * 2
    ampy = abs(ampy)
    ampy_neg = ampy * -1
    ms_dist2 = 0
    RR_list2 = []
    cnt = 0
    while (cnt < (len(peaklist2)-1)):
        # Calculate distance between beats in # of samples
        RR_interval2 = (peaklist2[cnt+1] - peaklist2[cnt])
        # Convert sample distances to ms distances
        ms_dist2 = ((RR_interval2 / fps) * 1000.0)
        RR_list2.append(ms_dist2)  # Append to list
        cnt += 1

    amp_medio1 = [ampx]
    amp_medio2 = [ampy]

    amp_medio1.extend(amp_medio2)

    amp_final = statistics.mean(amp_medio1)

    if not any((RR_list, RR_list)):
        print('não lista')
        bps = 0
        ms_dist = 'indetectable'
    else:
        # print(ms_dist)
        # print(ms_dist2)

        ms_dist = [ms_dist]
        ms_dist.append(ms_dist2)
        ms_dist = statistics.mean(ms_dist)
        # print(ms_dist)

        RR_list.extend(RR_list2)

        # 60000 ms (1 minute) / average R-R interval of signal
        bps = 60000 / np.mean(RR_list) / 60

    picos = (len(RR_list)/2) + 1

    if ms_dist != 'indetectable':
        value = ms_dist/1000
        value2 = str(round(value, 2))
    else:
        value2 = 'indetectable'

    # frequancia da batida
    duration = totaltime / picos
    value22 = 1/duration

    texto1 = "Beating frequency [Hz]: " + str(round(value22, 2))
    texto2 = 'Average of tail´s movement amplitude [degrees]: ' + str(
        round(amp_final, 2))
    texto2_1 = 'Sinoids detected: ' + str(picos)
    texto3 = 'FPS detected: ' + str(fps)

    text_6 = 'Beating cicle duration [seconds]: ' + str(round(
        duration, 2)) + '\n*Peaks with no spots were considered outliers and no computed.'

    plt.title(skip_state + '\n\n Detected angle x frame number')
    plt.xlim(0, total_frames)
    plt.plot(dataset, alpha=0.5, color='blue',
             label="raw signal")  # Plot semi-transparent HR
    plt.plot(mov_avg, color='green')  # Plot moving average
    # Plot detected peaks positive   %bps
    plt.scatter(peaklist1, ybeat1, color='red', label='_nolegend_')
    # Plot detected peaks negative
    plt.scatter(peaklist2, ybeat2, color='orange', label='_nolegend_')

    plt.legend(loc=4, framealpha=0.6)

    plt.xlabel('Absolute number of video frames \n\n Results:\n' + texto1 +
               '\n' + texto2 + '\n' + texto2_1 + '\n' + text_6 + '\n' + texto3 + '\n' + 'Coeficient of variation: positive is ' + str(round(cv_ybeat1, 2)) + ', negative is ' + str(round(cv_ybeat2, 2)) + ' and general is ' + str(round(cv_ybeat3, 2)))

    plt.ylabel('Tail angle (degrees)')

    plt.savefig(out_dir + 'line_plot.pdf', bbox_inches="tight", pad_inches=1,
                transparent=False, edgecolor='w', orientation='landscape')
    plt.savefig(out_dir + 'line_plot.png', bbox_inches="tight", pad_inches=1,
                transparent=False, edgecolor='w', orientation='landscape')

    # plt.show()
lista_original_filter = pd.read_csv(in_dir)

lista_original_filter = lista_original_filter.iloc[start_frame:finish_frame]

lista_original_filter['angle'].loc[(
    lista_original_filter['angle'].isnull())] = 0
lista_original_filter['angle'].loc[(lista_original_filter['angle'] > 40)] = 0
lista_original_filter['angle'].loc[(lista_original_filter['angle'] < -40)] = 0
no_zeros_ini = lista_original_filter[lista_original_filter != 0.00]
no_zeros = no_zeros_ini.dropna()
no_zeros.reset_index(drop=True, inplace=True)


no_zeros_list = no_zeros['angle'].tolist()
final_list = []
for idx, n in enumerate(no_zeros_list):
    if idx == 0:
        final_list.append(n)
    else:
        if n == no_zeros_list[idx - 1]:
            continue
        else:
            final_list.append(n)


peaks(final_list, fps, out_dir)
