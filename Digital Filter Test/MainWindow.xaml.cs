using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Digital_Filter_Test
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, RoutedEventArgs e)
        {
            chartCanvas.Children.Clear();
            SimpleSinusoidalTest();
        }


        private void SimpleSinusoidalTest()
        {
            Polyline pl;
            // Draw sine curve:
            pl = new Polyline();
            pl.Stroke = Brushes.Black;
            for (int i = 0; i < 360; i++)
            {
                double x = 2.0 * Math.PI * ((double)i / 360.0);
                double y = Math.Sin(x);
                pl.Points.Add(CurvePoint(
                new Point(x, y)));
            }
            chartCanvas.Children.Add(pl);
            // Draw cosine curve:
            pl = new Polyline();
            pl.Stroke = Brushes.Black;
            pl.StrokeDashArray = new DoubleCollection(
            new double[] { 4, 3 });

            for (int i = 0; i < 360; i++)
            {
                double x = 2.0 * Math.PI * ((double)i / 360.0);
                double y = Math.Cos(x);
                pl.Points.Add(CurvePoint(
                new Point(x, y)));
            }
            chartCanvas.Children.Add(pl);
        }

        private Point CurvePoint(Point pt)
        {
            const double xmin = 0;
            const double xmax = 2.0 * Math.PI;
            //
            const double ymin = -1.1;
            const double ymax = 1.1;
            //
            Point result = new Point();
            result.X = pt.X / xmax * chartCanvas.Width;
            result.Y = chartCanvas.Height - (pt.Y - ymin) / (ymax - ymin) * chartCanvas.Height;
            return result;
        }



        Random rand = new Random();
        Pen drawingPen = new Pen(Brushes.Black, 1);

        protected Point[] points;



        private void button2_Click(object sender, RoutedEventArgs e)
        {
            List<double> vals = new List<double>();


            // f(x) = (x3 + 3x2 − 6x − 8)/4
            for (int x = -100; x < 100; x++)
            {
                //double f = (Math.Pow(x, 3) + 3 * Math.Pow(x, 2) - 6 * x - 8) / 4;
                double f = Math.Pow(x, 2) + 10;
                vals.Add(f);

            }

            DrawFunction(vals);

        }

        private void DrawFunction(List<double> dvals)
        {
            double hmin = double.MaxValue;
            double hmax = double.MinValue;

            for (int x = 0; x < dvals.Count; x++)
            {
                double f = dvals[x];

                hmin = hmin > f ? hmin = f : hmin;
                hmax = hmax < f ? hmax = f : hmax;
            }


            double deltah = hmax - hmin;
            double stepsizeh = chartCanvas.ActualHeight / deltah;

            points = new Point[dvals.Count];

            double w = chartCanvas.ActualWidth / dvals.Count;
            double h = 1.0;

            double h_offset = chartCanvas.ActualHeight - stepsizeh * Math.Abs(hmin);


            points = new Point[dvals.Count];
            for (int i = 0; i < dvals.Count; i++)
            {
                points[i].X = (double)i;
                points[i].Y = h_offset - dvals[i] * stepsizeh * 0.9; // scale with 0.9 to leave borders at up&bottom.
            }



            chartCanvas.Children.Clear();

            for (int i = 1; i < dvals.Count; i++)
            {

                Line newLine = new Line();
                newLine.Stroke = Brushes.Black;
                newLine.Fill = Brushes.Black;
                newLine.StrokeLineJoin = PenLineJoin.Bevel;

                newLine.X1 = points[i - 1].X * w;
                newLine.Y1 = points[i - 1].Y;

                newLine.X2 = points[i].X * w;
                newLine.Y2 = points[i].Y;

                newLine.StrokeThickness = 1;
                chartCanvas.Children.Add(newLine);

            }


        }


        // SIGNAL
        List<double> v = new List<double>();


        private void button3_Click(object sender, RoutedEventArgs e)
        {
            v = new List<double>();

            // f(x) = (x3 + 3x2 − 6x − 8)/4


            /*
            Y(t) = (4/pi) sin (2 pi 100t) + (4/(3pi) ) * sin(2 pi 300 t) + (4/(5pi) ) * sin(2 pi 500 t)

            spectrum:
            100Hz  =>  1.3  volts signal amplitude
            300Hz  =>  0.45 volts signal amplitude
            500Hz  =>  0.30 volts signal amplitude

            Sampling Time delta = 10 micro seconds
            Total Sampling Time = 20 milli seconds
            Sampling Count = 2000

            Sampling Frequency = 1 MHz

            */


            double t = 0;
            double ts = 1 / 1.250E6;


            for (int x = 0; x < 5000; x++)
            {
                t = x * ts;

                double angle1 = 2.0 * Math.PI * 1000 * t;
                double angle2 = 2.0 * Math.PI * 20000 * t;
                double angle3 = 2.0 * Math.PI * 30000 * t;

                double y =
                      (4.0 / Math.PI) * Math.Sin(angle1)
                    + (4.0 / (3.0 * Math.PI)) * Math.Sin(angle2)
                    + (4.0 / (5.0 * Math.PI)) * Math.Sin(angle3);

                v.Add(y);
            }

            DrawFunction(v);

        }

        private void button4_Click(object sender, RoutedEventArgs e)
        {
            if (v.Count == 0)
            {
                MessageBox.Show("Create Signal by clicking 'SIGNAL' button");
            }

            double fc = 15E3;   // cutoff freq
            double TS = 1.25E6;  // Sampling rate


            double omega = 2 * Math.PI * (fc) * (1 / TS);
            double cos_omega = Math.Cos(omega);

            double alfa = (2 - Math.Sqrt(4 - 4 * cos_omega * cos_omega)) / (2.0 * cos_omega);
            double alfa2 = (2 + Math.Sqrt(4 - 4 * cos_omega * cos_omega)) / (2.0 * cos_omega);

            double f1 = (1 - alfa) / 2;
            double f2 = alfa;

            List<double> newvals = new List<double>();

            double prev_y = 0;
            double prev_x = v.Count > 0 ? v[0] : 0;

            for (int x = 1; x < v.Count; x++)
            {

                double y = f1 * v[x] + f1 * prev_x + f2 * prev_y;

                newvals.Add(y);

                prev_x = v[x];
                prev_y = y;
            }

            DrawFunction(newvals);

        }



        private static double[] FIR(double[] b, double[] x)
        {
            int M = b.Length;
            int n = x.Length;
            //y[n]=b0x[n]+b1x[n-1]+....bmx[n-M]
            var y = new double[n];
            for (int yi = 0; yi < n; yi++)
            {
                double t = 0.0;
                for (int bi = M - 1; bi >= 0; bi--)
                {
                    if (yi - bi < 0) continue;

                    t += b[bi] * x[yi - bi];
                }
                y[yi] = t;
            }
            return y;
        }


        double[] fir_coefs = 
        {
  2.2694038241719277,
  -8.201823869911099,
  11.908796635455738,
  -8.201823869911099,
  2.2694038241719277

        };


        // FIR
        // http://t-filter.engineerjs.com/

        private void button5_Click(object sender, RoutedEventArgs e)
        {
            if (v.Count == 0)
            {
                MessageBox.Show("Create Signal by clicking 'SIGNAL' button");
            }


            List<double> newvals = new List<double>();


            double[] x = v.ToArray();

            // FIR(double[] b, double[] x)
            double[] k = FIR(fir_coefs, x);

            for (int n = 0; n < k.Length; n++)
                newvals.Add(k[n]);


            DrawFunction(newvals);


        }


    }
}
