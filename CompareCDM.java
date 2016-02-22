/* CompareCDM.java */

/* This file is used to compare the coordinate descend method with CG method
 * examples here using coordinate descend method
 */
 

import java.io.*;
import java.util.*;

public class CompareCDM {
    
    public static void main(String[] args) {
        
        // Open a new file "CDM.txt", write the title into the first line of file.
        try {
            FileWriter fw = new FileWriter("CDM.txt");
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("x0\t\t\t y0\t\t\t final x \t\t final y\t\t final f(x,y) \t timesOfEvalution \t timer(ns)");
            bw.newLine();
            bw.close();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println("error" + e.getMessage());
        }
        
    	/* 
    	 * random generate a hundred set of (x0,y0) using seed 10062015(due date),
    	 * and run a hundred times optimization.
    	 */
    	Random random = new Random (10062015);
    	int sampleNumber = 100; // user define sample number.
    	double[] x0 = new double [sampleNumber];
    	double[] y0 = new double [sampleNumber];
    	double[] xf = new double [sampleNumber];
    	double[] yf = new double [sampleNumber];
    	double[] ff = new double [sampleNumber];
    	double[] toe = new double [sampleNumber];
    	double[] timer = new double [sampleNumber];

    	for (int i = 0; i <= sampleNumber-1; i++){

    		// timer starts
    		long tStart = System.nanoTime();

    		// starting point x0 is randomly chosen from 0 ~ +50 (user choose)
    		double xStart = random.nextDouble() *1000-500;
    		double yStart = random.nextDouble() *1000-500;

    		Optimization opt = new Optimization(xStart, yStart);

    		double [] resultX = new double[3];
    		double [] resultY = new double[3];
    		double tmpX, tmpY, dis1, dis2;
    		toe[i] = 0;
    		do {
    			tmpX = opt.xCoordinate;
            	tmpY = opt.yCoordinate;

            	resultX = opt.quadFitX(); // apply quadFitX
            	opt.xCoordinate = resultX[0];
                toe[i]++;

            	resultY = opt.quadFitY(); // apply quadFitY
            	opt.yCoordinate = resultY[0];
            	toe[i]++;

            	dis1 = Math.pow((Math.pow(opt.xCoordinate,2)+Math.pow(opt.yCoordinate,2)), 0.5);
            	dis2 = Math.pow((Math.pow(tmpX,2)+Math.pow(tmpY,2)), 0.5);

            } while (Math.abs(dis1 - dis2) > 0.000001);

            // store the x0, final x, final f(x), and times of evaluation.
            x0[i] = xStart;
            y0[i] = yStart;
            xf[i] = resultX[0];
            yf[i] = resultY[0];
            ff[i] = resultY[1];

            // timer ends
            double tEnd = System.nanoTime();
            double tDelta = tEnd - tStart;
            timer[i] = tDelta;

            // write the results into file.
            try {
                FileWriter fw = new FileWriter("2dSampleAnalysis.txt",true);
                BufferedWriter bw = new BufferedWriter(fw);
                bw.write(xStart+"\t"+yStart+"\t"+xf[i]+"\t"+yf[i]+"\t"+ff[i]+"\t"+toe[i]+"\t"+tDelta);
                bw.newLine();
                bw.close();
            } catch(Exception e)
            {
                e.printStackTrace();
                System.out.println("error"+e.getMessage());
            }
        }
        
        System.out.printf("Open 2dSampleAnalysis.txt in current directory to see.\n");
        
        /*
         * Mathematical analysis including average, variance and error bar for x, y, f(x), times of evolution,
         * running time, and relative distance.
         */
        // average, variance, and error bar of x
        double sumXf = 0;
        for (int j = 0; j<=xf.length - 1; j++) {
    	   sumXf = sumXf + xf[j];
        }
        double avrXf = sumXf/xf.length;
        System.out.printf("the average x of 100 samples is %.15g\n", avrXf);
        
        double varXf = 0;
        for (int j = 0; j<xf.length; j++) {
        	varXf = varXf+ Math.pow((xf[j] - avrXf), 2);
        }
        varXf = varXf/xf.length;
        System.out.printf("the sample variance of x is %.15g\n", varXf);
        
        double errorXf = Math.pow(varXf/(sampleNumber-1), 0.5);
        System.out.printf("the sample error bar of x is %.15g\n\n", errorXf);
        
        // average, variance and error bar of y

        double sumYf = 0;
        for (int j = 0; j<=yf.length - 1; j++) {
    	   sumYf = sumYf + yf[j];
        }
        double avrYf = sumYf/yf.length;
        System.out.printf("the average y of 100 samples is %.15g\n", avrYf);
        
        double varYf = 0;
        for (int j = 0; j<yf.length; j++) {
        	varYf = varYf+ Math.pow((yf[j] - avrYf), 2);
        }
        varYf = varYf/yf.length;
        System.out.printf("the sample variance of y is %.15g\n", varYf);
        
        double errorYf = Math.pow(varYf/(sampleNumber-1), 0.5);
        System.out.printf("the sample error bar of y is %.15g\n\n", errorYf);
        
        
        // average, variance, and error bar of f(x) 
        double sumF = 0;
        for (int j = 0; j<ff.length; j++) {
        	sumF = sumF + ff[j];
        }
       double avrF = sumF/ff.length;
       System.out.printf("the average of f of 100 samples is %.15g\n", avrF);
       
       double varF = 0;
       for (int j = 0; j<ff.length; j++) {
       	varF = varF+ Math.pow((ff[j] - avrF), 2);
       }
       varF = varF/ff.length;
       System.out.printf("the sample variance of f(x,y) is %.15g\n", varF);  
       
       double errorF = Math.pow(varF/(sampleNumber-1), 0.5);
       System.out.printf("the sample error bar of f(x,y) is %.15g\n\n", errorF);
       
       
       // average, variance, and error bar of times of evaluation
       double sumToe = 0;
       for (int j = 0; j<toe.length; j++) {
    	   sumToe = sumToe + toe[j];
       }
       double avrToe = sumToe/toe.length;
       System.out.printf("the average of times of evaluation is %.0f\n ", avrToe);
       
       double varToe = 0;
       for (int j = 0; j<toe.length; j++) {
       	varToe = varToe+ Math.pow((toe[j] - avrToe), 2);
       }
       varToe = varToe/toe.length;
       System.out.printf("the sample variance of times of evaluation is %.15g\n", varToe);  
       
       double errorToe = Math.pow(varToe/(sampleNumber-1), 0.5);
       System.out.printf("the sample error bar of times of evaluation is %.15g\n\n", errorToe);
       
       
       // average of wall clock running time
       double sumTimer = 0;
       for (int j= 0; j<timer.length; j++) {
    	   sumTimer = sumTimer + timer[j];
       }
       double avrTimer = sumTimer/timer.length;
       System.out.printf("the average of running time is %f nano seconds\n", avrTimer);

       double varTimer = 0;
       for (int j = 0; j<timer.length; j++) {
       	varTimer = varTimer+ Math.pow((timer[j] - avrTimer), 2);
       }
       varTimer = varTimer/timer.length;
       System.out.printf("the sample variance of running time is %.15g\n", varTimer);  
       
       double errorTimer = Math.pow(varTimer/(sampleNumber-1), 0.5);
       System.out.printf("the sample error bar of running time is %.15g\n\n", errorTimer);
       
       
       
       // average of relative distance between output optimum and true mathematical optimum
       /*
        * f = 1/2 * (x1^2 + x2^2) + 1/2 * (x1+x2-1)^2 - (x1 + 2*x2);
        * mathematical optimum is f = -11/6  at x1=1/3, x2=4/3.
        */
       double optimum = -11/6;
       double sumDis = 0;
       for (int j= 0; j<ff.length; j++) {
    	   sumDis = sumDis + Math.abs(ff[j] - optimum);
       }
       double avrDis = sumDis/ff.length;
       System.out.printf("average of relative distance is %.15g\n", avrDis);

       double varDis = 0;
       for (int j = 0; j<ff.length; j++) {
       	varDis = varDis+ Math.pow((Math.abs(ff[j] - optimum) - avrDis), 2);
       }
       varDis = varDis/ff.length;
       System.out.printf("the sample variance of relative distance is %.15g\n", varDis);  
       
       double errorDis = Math.pow(varDis/(sampleNumber-1), 0.5);
       System.out.printf("the sample error bar of relative distance is %.15g\n\n", errorDis);

       
       // write the results into file.
       try {
           FileWriter fw = new FileWriter("CDM.txt",true);
           BufferedWriter bw = new BufferedWriter(fw);
           bw.write("\n"+avrXf+"\t"+varXf+"\t"+errorXf+"\n"+avrYf+"\t"+varYf+"\t"+errorYf+"\n"
                    +avrF+"\t"+varF+"\t"+errorF+"\n"+avrToe+"\t"+varToe+"\t"+errorToe+"\n"
                    +avrTimer+"\t"+varTimer+"\t"+errorTimer+"\n"+avrDis+"\t"+varDis+"\t"+errorDis+"\t");
           bw.newLine();
           bw.close();
       } catch(Exception e)
       {
           e.printStackTrace();
           System.out.println("error"+e.getMessage());
       }
    }
}