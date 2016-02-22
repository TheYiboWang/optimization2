/* Optimization.java */

/* The Optimization class is the process of 1d optimization.
 *
 * First use bracketing to determine x1,x2,x3(or y1, y2, y3) where f1>f2<f3.
 * Second use quadFitX or quadFitY to determine the x(or y) where f is minimum.
 * Stopping criteria is determined by user in this program.
 * The saftguard and golden section step are used to make sure quadFit works well
 *
 * @author Yibo Wang, any questions please contact yibow6@uci.edu
 * 10-06-2015
 */
//import Function.class;

public class Optimization {
    double xCoordinate;
    double yCoordinate;
    double zCoordinate;
    Function fun;
    
    public Optimization(double x0) {
        xCoordinate = x0;
        yCoordinate = 0;
        zCoordinate = 0;
        fun = new Function(x0);
    }
    
    public Optimization(double x0, double y0) {
        xCoordinate = x0;
        yCoordinate = y0;
        zCoordinate = 0;
        fun = new Function(x0, y0);
    }
    
    public Optimization(double x0, double y0, double z0) {
        xCoordinate = x0;
        yCoordinate = y0;
        zCoordinate = z0;
        fun = new Function(x0, y0, z0);
    }
    
    /*
     * bracketing and quadratic fit method on x axis
     */
    public double[] quadFitX() {
        
        // apply three point pattern bracketing to get f1>f2<f3 //
        double delta = 0.000001; // this value is self-defined.
        double TAU = 1.618; // book suggested.
        
        double[] x = new double [5]; // an array stores x1, x2, x3, x4. the 0 position is not used.
        double[] f = new double [5]; // an array stores f1, f2, f3, f4. the 0 position is not used.
        
        x[1] = fun.xCoordinate;
        f[1] = fun.calculateF(); // calculate f1 and store
        x[2] = x[1] + delta; // calculate x2 and store
        fun.xCoordinate = x[2]; // calculate f2 and store
        f[2] = fun.calculateF();
        
        // determine f2 and x2
        if (f[2] > f[1]) {
            //rename and swapping
            double tmp = f[1]; f[1] = f[2]; f[2] = tmp;
            tmp = x[1]; x[1] = x[2]; x[2] = tmp;
            delta = -delta;
        }
        
        // determine f3 and x3
        delta = TAU * delta;
        x[3] = x[2] + delta;
        fun.xCoordinate = x[3];
        f[3] = fun.calculateF();
        
        // update f3 till f2<f3
        while (f[2] > f[3]) {
            //rename f2 as f1, f3 as f2, x2 as x1, x3 as x2
            f[1] = f[2]; f[2] = f[3];
            x[1] = x[2]; x[2] = x[3];
            
            delta = TAU * delta;
            x[3] = x[2] + delta;
            fun.xCoordinate = x[3];
            f[3] = fun.calculateF();
        }
        
        // set x1<x2<x3
        if (x[1] > x[3]) {
            double tmp = x[1]; x[1] = x[3]; x[3] = tmp;
            tmp = f[1]; f[1] = f[3]; f[3] = tmp;
        }
        
        
        //  quadratic fit method for x axis
        double xtol = 0.000001; // stopping criteria, user defined.
        double countf4 = 0.0; // times of evaluation.
        
        // calculate x4 and f4
        while( Math.abs(x[3] - x[1]) > 2.0 * xtol) {
            
            double a = (x[1] - x[2]) * (x[1] - x[3]);
            double b = (x[2] - x[1]) * (x[2] - x[3]);
            double c = (x[3] - x[1]) * (x[3] - x[2]);
            x[4] = 0.5 * (f[1] * (x[2] + x[3])/a + f[2] * (x[1] + x[3])/b + f[3] * (x[1] + x[2])/c) / (f[1]/a + f[2]/b + f[3]/c);
            
            // safeguard x4 cannnot be too close to any of the x1, x2, and x3.
            double measure = 1/8 * Math.min(Math.abs(x[3]-x[2]), Math.abs(x[2]-x[1]));
            if (Math.abs(x[4]-x[1]) < measure) {
                x[4] = x[1] + measure;
            }
            if (Math.abs(x[4]-x[3]) < measure) {
                x[4] = x[3] - measure;
            }
            if (Math.abs(x[4]-x[2]) < measure && x[2] > 0.5 * (x[1]+x[3])) {
                x[4] = x[2] - measure;
            }
            if (Math.abs(x[4]-x[2]) < measure && x[2] <= 0.5 * (x[1]+x[3])) {
                x[4] = x[2] + measure;
            }
            
            // golden section step is impletemented if reduce ratio is small
            if (x[2] <= (x[1] + x[3])/2) {
                x[4] = x[2] + (1-0.618) * (x[3] - x[2]);
            } else {
                x[4] = x[3] - (1-0.618) * (x[2] - x[1]);
            }
            
            // calculate f(x) and number of evalutions increase by 1.
            fun.xCoordinate = x[4];
            f[4] = fun.calculateF();
            countf4++;
            
            // renew x1, x2, and x3
            if (x[4] > x[2] && f[4] >= f[2]) {
                x[3] = x[4];
                f[3] = f[4];
            }
            if (x[4] > x[2] && f[4] < f[2]) {
                x[1] = x[2]; x[2] = x[4];
                f[1] = f[2]; f[2] = f[4];
            }
            if (x[4] < x[2] && f[4] >= f[2]) {
                x[1] = x[4];
                f[1] = f[4];
            }
            if (x[4] < x[2] && f[4] < f[2]) {
                x[3] = x[2]; x[2] = x[4];
                f[3] = f[2]; f[2] = f[4];
            }
        }
        
        return new double[] {x[4], f[4], countf4};
    }
    
    /*
     * bracketing and quadratic fit method on y axis
     */
    public double[] quadFitY() {
        
        /* apply three point pattern bracketing to get f1>f2<f3 */
        double delta = 0.1; // this value is user-defined.
        double TAU = 1.618; // book suggested.
        
        double[] y = new double [5]; // an array stores y1, y2, y3, y4. the 0 position is not used.
        double[] f = new double [5]; // an array stores f1, f2, f3, f4. the 0 position is not used.
        
        y[1] = fun.yCoordinate;
        f[1] = fun.calculateF(); // calculate f1 and store
        y[2] = y[1] + delta; // calculate x2 and store
        fun.yCoordinate = y[2]; // calculate f2 and store
        f[2] = fun.calculateF();
        
        // determine f2 and y2
        if (f[2] > f[1]) {
            //rename and swapping
            double tmp = f[1]; f[1] = f[2]; f[2] = tmp;
            tmp = y[1]; y[1] = y[2]; y[2] = tmp;
            delta = -delta;
        }
        
        // determine f3 and y3
        delta = TAU * delta;
        y[3] = y[2] + delta;
        fun.yCoordinate = y[3];
        f[3] = fun.calculateF();
        
        // update f3 till f2<f3
        while (f[2] > f[3]) {
            //rename f2 as f1, f3 as f2, y2 as y1, y3 as y2
            f[1] = f[2]; f[2] = f[3];
            y[1] = y[2]; y[2] = y[3];
            
            delta = TAU * delta;
            y[3] = y[2] + delta;
            fun.yCoordinate = y[3];
            f[3] = fun.calculateF();
        }
        
        // set y1<y2<y3
        if (y[1] > y[3]) {
            double tmp = y[1]; y[1] = y[3]; y[3] = tmp;
            tmp = f[1]; f[1] = f[3]; f[3] = tmp;
        }
        
        
        //  quadratic fit method for y axis
        double ytol = 0.000001; // stopping criteria, user defined.
        double countfy4 = 0.0; // times of evaluation.
        
        // calculate y4 and f4
        while( Math.abs(y[3] - y[1]) > 2.0 * ytol) {
            
            double a = (y[1] - y[2]) * (y[1] - y[3]);
            double b = (y[2] - y[1]) * (y[2] - y[3]);
            double c = (y[3] - y[1]) * (y[3] - y[2]);
            y[4] = 0.5 * (f[1] * (y[2] + y[3])/a + f[2] * (y[1] + y[3])/b + f[3] * (y[1] + y[2])/c) / (f[1]/a + f[2]/b + f[3]/c);
            
            // safeguard x4 cannnot be too close to any of the y1, y2, and y3.
            double measure = 1/8 * Math.min(Math.abs(y[3]-y[2]), Math.abs(y[2]-y[1]));
            if (Math.abs(y[4]-y[1]) < measure) {
                y[4] = y[1] + measure;
            }
            if (Math.abs(y[4]-y[3]) < measure) {
                y[4] = y[3] - measure;
            }
            if (Math.abs(y[4]-y[2]) < measure && y[2] > 0.5 * (y[1]+y[3])) {
                y[4] = y[2] - measure;
            }
            if (Math.abs(y[4]-y[2]) < measure && y[2] <= 0.5 * (y[1]+y[3])) {
                y[4] = y[2] + measure;
            }
            
            // golden section step is impletemented if reduce ratio is small
            if (y[2] <= (y[1] + y[3])/2) {
                y[4] = y[2] + (1-0.618) * (y[3] - y[2]);
            } else {
                y[4] = y[3] - (1-0.618) * (y[2] - y[1]);
            }
            
            // calculate f(y) and number of evalutions increase by 1.
            fun.yCoordinate = y[4];
            f[4] = fun.calculateF();
            countfy4++;
            
            // renew y1, y2, and y3
            if (y[4] > y[2] && f[4] >= f[2]) {
                y[3] = y[4];
                f[3] = f[4];
            }
            if (y[4] > y[2] && f[4] < f[2]) {
                y[1] = y[2]; y[2] = y[4];
                f[1] = f[2]; f[2] = f[4];
            }
            if (y[4] < y[2] && f[4] >= f[2]) {
                y[1] = y[4];
                f[1] = f[4];
            }
            if (y[4] < y[2] && f[4] < f[2]) {
                y[3] = y[2]; y[2] = y[4];
                f[3] = f[2]; f[2] = f[4];
            }
        }
        
        return new double[] {y[4], f[4], countfy4};
    }
}