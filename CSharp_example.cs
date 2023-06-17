using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.NetworkInformation;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using System.Threading.Tasks;
using static Csharp_test.Program;

namespace Csharp_test
{
    public class OWON
    {
        public enum wavelets
        {
            hermitian1 = 0,
            hermitian2 = 1,
            hermitian3 = 2,
            hermitian4 = 3,
            poisson2 = 4,
            morlet = 5,
            modified_morlet = 6
        };

        public enum window
        {
            none,
            square,
            triangle,
            RCF // Raised Cosine Filter
        };

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        public static extern int getOWONData2(
            int sample_size,
            double sample_step,
            int sample_perSec_max,
            double voltage_max,
            int trigger_level,
            int reading_number,
            double[] result_arr,
            char[] result_fileName
        );

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        unsafe public static extern void sizing_toDouble(
            int t_input_startIndex,
            int t_input_endIndex,
            double t_input_step,
            int t_input_size,
            int t_sizing,
            double* t_sizing_start,
            double* t_sizing_end,
            double* t_sizing_step,
            int* t_sizing_size
        );

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        public static extern void make_waveletFunction_equalStep(double t_start, double t_step, int t_size, double f_start, double f_step, int f_size, wavelets wavelet, double[,] res);

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        public static extern void waveletTransform_fromIndexToIndex(
            char[] inputFIleName,
            double t_input_step, double t_input_size,
            int t_start_index, int t_end_index, int t_sizing,
            double f_start, double f_step, int f_size,
            double[,] psitaus,
            string outputFile_name,
            OWON.window win,
            double[] win_top_t,
            double[] win_top_f,
            int win_top_n,
            double[] win_bot_t,
            double[] win_bot_f,
            int win_bot_n
            );

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        public static extern void waveletTransform_fromTimeToTime(
            char[] inputFileName,
            double t_input_step, double t_input_size,
            int t_start, int t_end, int t_sizing,
            double f_start, double f_step, int f_size,
            double[,] psitaus,
            char[] outputFile_name,
            window win,
            double[] win_top_t,
            double[] win_top_f,
            int win_top_n,
            double[] win_bot_t,
            double[] win_bot_f,
            int win_bot_n
        );

        [DllImport("C:/Users/user1/source/repos/owon_signalProcessing/x64/Release/owon_signalProcessing.dll")]
        public static extern void backWavelet(char[] transformedSignal_fileName,
            double t_start, double t_step, int t_size,
            double f_start, double f_step, int f_size,
            double[,] waveletFunction,
            char[] originalSignal_fileName
            );
    }

    class Program
    {
        static void Main(string[] args)
        {
            //Variables for getting data from OWON
            int sample_size = 1000;
            double sample_step = 10E-9; //0 to min
            int sample_perSec_max = 125000000; //0 to default (125E6)
            double voltage_max = 10;
            int reading_number = 10;
            int trigger_level = 0; //from -5 to 5, where 5 equals voltage_max
            double[] read_arr = new double[sample_size];
            string fileName = "read.txt";

            fileName += '\0';
            char[] read_fileName = new char[fileName.Length];
            read_fileName = fileName.ToCharArray();

            for (int i = 0; i < sample_size; i++)
                read_arr[i] = 0;

            OWON.getOWONData2(sample_size, sample_step, sample_perSec_max, voltage_max, trigger_level, reading_number, read_arr, read_fileName);

            //Variables for wavelet transform
            int t_sizing = 1; //If you need to use only every t_sizing element of samples

            //Frequency splitting
            int f_start = 0;
            int f_step = 1000;
            int f_size = 200;

            //The wavelet that will be used in the transform
            OWON.wavelets wavelet = OWON.wavelets.modified_morlet;

            OWON.window win = OWON.window.RCF;
            int win_top_n = 0;
            double[] win_top_t = new double[win_top_n];
            double[] win_top_f = new double[win_top_n];

            int win_bot_n = 2;
            double[] win_bot_t = new double[win_bot_n];
            double[] win_bot_f = new double[win_bot_n];
            win_bot_t[0] = 0 * sample_step;
            win_bot_t[1] = sample_size * sample_step;

            win_bot_f[0] = 20 * f_step;
            win_bot_f[1] = 20 * f_step;

            //Let it be "" if you don't need to transform the signal
            //string transformedSignal_fileName = "";
            string transformedSignal_fileName = "signal_wavelet.txt";
            string backTransformedSignal_fileName = "signal_backWavelet.txt";

            transformedSignal_fileName += '\0';
            char[] char_transformedSignal_fileName = new char[transformedSignal_fileName.Length];
            char_transformedSignal_fileName = transformedSignal_fileName.ToCharArray();

            backTransformedSignal_fileName += '\0';
            char[] char_backTransformedSignal_fileName = new char[backTransformedSignal_fileName.Length];
            char_backTransformedSignal_fileName = backTransformedSignal_fileName.ToCharArray();

            if (transformedSignal_fileName != "")
            {
                double t_start = 0, t_end = 0, t_step = 0;
                int t_size = 0;

                //It is used in backWavelet() and in waveletTransform_fromTimeToTime()
                unsafe {
                    OWON.sizing_toDouble(0, sample_size, sample_step, sample_size, t_sizing,
                        &t_start, &t_end, &t_step, &t_size
                        );
                }

                //Wavelet function generation
                double[,] wavelet_function = new double[f_size, (2 * t_size)];
                OWON.make_waveletFunction_equalStep(t_start, t_step, t_size, f_start, f_step, f_size, wavelet, wavelet_function);

                //Wavelet transforms
                //double[,] wavelet_transform = new double[f_size, t_size];
                OWON.waveletTransform_fromIndexToIndex(read_fileName, sample_step, sample_size, 0, sample_size, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName, win, win_top_t, win_top_f, win_top_n, win_bot_t, win_bot_f, win_bot_n);
                //wavelet_trans::waveletTransform_fromTimeToTime(read_fileName, sample_step, sample_size, t_start, t_end, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName);
                //printf("The signal has been transformed.\n");

                if (backTransformedSignal_fileName != "")
                    OWON.backWavelet(char_transformedSignal_fileName, t_start, t_step, t_size, f_start, f_step, f_size, wavelet_function, char_backTransformedSignal_fileName);
                //printf("The wavelet has been retransformed.\n");
            }
        }
    }
}
