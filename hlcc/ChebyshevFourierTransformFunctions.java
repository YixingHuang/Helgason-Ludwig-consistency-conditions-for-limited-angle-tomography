package edu.stanford.rsl.tutorial.hlcc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix.MatrixNormType;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import WeightedTV.wtv2D.FanBeamProjector2D;
import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix.InversionType;
import edu.stanford.rsl.tutorial.fan.CosineFilter;
import WeightedTV.wtv2D.FanBeamBackprojector2D;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import edu.stanford.rsl.tutorial.interpolation.FanParameters;
import edu.stanford.rsl.tutorial.phantoms.Phantom;
import edu.stanford.rsl.tutorial.phantoms.SheppLogan;
import edu.stanford.rsl.tutorial.wedgefilter.DoubleWedgeFilterFanES;
import ij.ImageJ;
import ij.IJ;
import ij.ImagePlus;

public class ChebyshevFourierTransformFunctions {
	public boolean extendToNextPowerOfTwo = false;
	
	public Grid2DComplex sinogramFourierError(Grid2DComplex sinoFourierIdeal, Grid2DComplex sinoFourier){
		Grid2DComplex sinoFourierError = sinoFourierIdeal.clone();
		sinoFourier.getGridOperator().subtractBy(sinoFourierError, sinoFourier);
		return sinoFourierError;
		}
	
	public Grid2DComplex FourierTransform(Grid2D sinogram){
		Grid2DComplex cCopy  = new Grid2DComplex(sinogram,false);		
		cCopy.transformForward();
		// FFT-shift the 2D-FT to fit the double wedge filter axis', where the zero frequency is centered
		cCopy.fftshift();
		return cCopy;
	}
	
	
	public double rmse(Grid2D recon, Grid2D recon_data) {
		double err = 0;
		Grid2D temp = new Grid2D(recon);
		NumericPointwiseOperators.subtractBy(temp, recon_data);
		err = Grid2Dnorm(temp);
		err = err / (temp.getSize()[0] * temp.getSize()[1]);
		err = Math.sqrt(err);
		return err;
	}
	

	private double Grid2Dnorm(Grid2D recon) {
		
		double d = 0;
		for (int row = 0; row < recon.getSize()[0]; row++)
			for (int col = 0; col < recon.getSize()[1]; col++)
				d = d + recon.getAtIndex(row, col) * recon.getAtIndex(row, col);
		return d;
	}
		
	/**
	 * Helgason-Ludwig Moment Transform
	 * @param sinogram
	 * @param k
	 * @return
	 */
	
	public Grid1D HelgasonLudwigMomentTransform(Grid2D sinogram, int k){
		Grid1D coefs;
		Grid1D momentSum = new Grid1D(sinogram.getSize()[1]);
		double deltaS = 2.f / sinogram.getSize()[0]; // here regard the detector from -1 to 1, maxS = 2, deltaS = 2 / detectorSize;
		
		Grid1D sPowerk = new Grid1D(sinogram.getSize()[0]);
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++)
			sPowerk.setAtIndex(i, (float) (-1 + (i + 0.5) * deltaS));  //FIXME
		
		sPowerk.getGridOperator().pow(sPowerk, k);
		
		Grid1D sub;
		for ( int i = 0; i < sinogram.getSize()[1]; i ++){
			sub = new Grid1D(sinogram.getSubGrid(i));
			sPowerk.getGridOperator().multiplyBy(sub, sPowerk);
			momentSum.setAtIndex(i, (float) sub.getGridOperator().sum(sub));
		}
			
		Grid1DComplex coeFourier = new Grid1DComplex (momentSum, false);
		coeFourier.transformForward();
		coefs = coeFourier.getMagSubGrid(0, coeFourier.getSize()[0]);
		
		
		return coefs;
	}
	
	/**
	 * Chebyshev Fourier transform coefficients
	 * @param sinogram
	 * @param n
	 * @return
	 */
	public Grid1DComplex ChebyshevFourierTransform(Grid2D sinogram, int n){
		Grid1D momentCurve = ChebyshevFourierTransformMoment(sinogram, n);
		Grid1DComplex coeFourier = new Grid1DComplex (momentCurve, extendToNextPowerOfTwo);
		coeFourier.transformForward();
		return coeFourier;
	}
	
	public Grid1DComplex ChebyshevFourierTransformFromMomentCurve(Grid1D momentCurve){
	
		Grid1DComplex coeFourier = new Grid1DComplex (momentCurve, extendToNextPowerOfTwo);
		coeFourier.transformForward();
		return coeFourier;
	}
	
	/**
	 * 
	 * 
	 * @param n the order of the moment curve
	 * @param numThetaLimited the number of limited angle projections, i. e., 
	 * the number of points measured on the moment curve
	 * @param numThetaPI the number of points in angular range PI
	 * @param momentCurve ground truth moment curve, 360 degree
	 * @return
	 */
	public Grid1D GPAlgorithmForMomentCurve(int n, int numThetaLimited, int numThetaPI, Grid1D momentCurve, int maxIter){
		Grid1D momentCurveLimited = new Grid1D(momentCurve);
		Grid1D momentCurveRestored;
		Grid1DComplex coeFourierLimited;
		
		for(int i = numThetaLimited; i < numThetaPI; i++){
			momentCurveLimited.setAtIndex(i, 0);
			momentCurveLimited.setAtIndex(i + numThetaPI, 0);
		}
	
		/*
		for(int i = 0; i < numThetaLimited; i++){
			momentCurveInitial.setAtIndex(i, momentCurve.getAtIndex(i));
			momentCurveInitial.setAtIndex(i + numThetaLimited, momentCurve.getAtIndex(i + numThetaPI));
		}
		coeFourierInitial = new Grid1DComplex (momentCurveInitial, extendToNextPowerOfTwo);
		coeFourierInitial.transformForward();
		consistencyConstrain(coeFourierInitial, n);
		coeFourierLimited = extendFourierLength(coeFourierInitial, momentCurve.getNumberOfElements());
		momentCurveLimited = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
		*/
		float err = 0, val;
		for(int l = 0; l < maxIter; l++){
			coeFourierLimited = new Grid1DComplex (momentCurveLimited, extendToNextPowerOfTwo);
			coeFourierLimited.transformForward();
			//consistencyConstrain(coeFourierLimited, n);
			consistencyConstrain(coeFourierLimited, n);
			momentCurveRestored = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
			err = 0;
			for(int i = numThetaLimited; i < numThetaPI; i++){
				momentCurveLimited.setAtIndex(i, momentCurveRestored.getAtIndex(i));
				val = momentCurveRestored.getAtIndex(i) - momentCurve.getAtIndex(i);
				err = err + val * val;
				momentCurveLimited.setAtIndex(i + numThetaPI, momentCurveRestored.getAtIndex(i+numThetaPI));
			}
			err = err/ (numThetaPI-numThetaLimited);
			err = (float)Math.sqrt(err);
			System.out.println("n = " + n + " error = " + err);
		}
		
	
		return momentCurveLimited;
	}
	
	public Grid1D GPAlgorithmForMomentCurveWithSoftThresholding(int n, int numThetaLimited, int numThetaPI, Grid1D momentCurve, int maxIter, float tau){
		Grid1D momentCurveLimited = new Grid1D(momentCurve);
		Grid1D momentCurveInitial = new Grid1D(2*numThetaLimited);
		Grid1D momentCurveRestored;
		Grid1DComplex coeFourierLimited, coeFourierInitial;
		
		for(int i = numThetaLimited; i < numThetaPI; i++){
			momentCurveLimited.setAtIndex(i, 0);
			momentCurveLimited.setAtIndex(i + numThetaPI, 0);
		}
	
		/*
		for(int i = 0; i < numThetaLimited; i++){
			momentCurveInitial.setAtIndex(i, momentCurve.getAtIndex(i));
			momentCurveInitial.setAtIndex(i + numThetaLimited, momentCurve.getAtIndex(i + numThetaPI));
		}
		coeFourierInitial = new Grid1DComplex (momentCurveInitial, extendToNextPowerOfTwo);
		coeFourierInitial.transformForward();
		consistencyConstrain(coeFourierInitial, n);
		coeFourierLimited = extendFourierLength(coeFourierInitial, momentCurve.getNumberOfElements());
		momentCurveLimited = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
		*/
		float err = 0, val;
		for(int l = 0; l < maxIter; l++){
			coeFourierLimited = new Grid1DComplex (momentCurveLimited, extendToNextPowerOfTwo);
			coeFourierLimited.transformForward();
			//consistencyConstrain(coeFourierLimited, n);
			consistencyConstrain(coeFourierLimited, n);
			softThresholdingOnFourierCoefficients(coeFourierLimited, tau);
			momentCurveRestored = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
			err = 0;
			for(int i = numThetaLimited; i < numThetaPI; i++){
				momentCurveLimited.setAtIndex(i, momentCurveRestored.getAtIndex(i));
				val = momentCurveRestored.getAtIndex(i) - momentCurve.getAtIndex(i);
				err = err + val * val;
				momentCurveLimited.setAtIndex(i + numThetaPI, momentCurveRestored.getAtIndex(i+numThetaPI));
			}
			err = err/ (numThetaPI-numThetaLimited);
			err = (float)Math.sqrt(err);
			System.out.println("n = " + n + " error = " + err);
		}
		
	
		return momentCurveLimited;
	}
	
	public void softThresholdingOnFourierCoefficients(Grid1DComplex coeFourier, float tau){
		for(int i = 0; i < coeFourier.getSize()[0]; i++){
			coeFourier.setImagAtIndex(i, softThresholding(coeFourier.getImagAtIndex(i), tau));
			coeFourier.setRealAtIndex(i, softThresholding(coeFourier.getRealAtIndex(i), tau));
		}
		
	}
	
	
	private float softThresholding(float val, float tau){
		if(val > tau)
			val = val - tau;
		else if (val < -tau)
			val = val + tau;
		else
			val = 0;
		return val;
	}
	private Grid1DComplex extendFourierLength(Grid1DComplex coeFourier, int len2){
		int len1 = coeFourier.getSize()[0];
		System.out.println("len1 = "+len1);
		Grid1DComplex coeFourier2 = new Grid1DComplex(len2);
		for(int i = 0; i <= len1/2; i++){
			coeFourier2.setImagAtIndex(i, coeFourier.getImagAtIndex(i));
			coeFourier2.setRealAtIndex(i, coeFourier.getRealAtIndex(i));
		}
		for(int i = 1; i <= len1; i++){
			coeFourier2.setImagAtIndex(len2-i, coeFourier.getImagAtIndex(len1-i));
			coeFourier2.setRealAtIndex(len2-i, coeFourier.getRealAtIndex(len1-i));
		}
		
		return coeFourier2;
		
	}
	/**
	 * 
	 * 
	 * @param n the order of the moment curve
	 * @param numThetaLimited the number of limited angle projections, i. e., 
	 * the number of points measured on the moment curve
	 * @param numThetaPI the number of points in angular range PI
	 * @param momentCurve ground truth moment curve, 360 degree
	 * @return
	 */
	public Grid1D GPAlgorithmForMomentCurve2(int n, int numThetaLimited, int numThetaPI, Grid1D momentCurve, int maxIter){
		Grid1D momentCurveLimited = new Grid1D(momentCurve);
		Grid1D momentCurveRestored, momentCurveRestored0;
		Grid1DComplex coeFourierLimited;
		for(int i = numThetaLimited; i < numThetaPI; i++){
			momentCurveLimited.setAtIndex(i, 0);
			momentCurveLimited.setAtIndex(i + numThetaPI, 0);
		}
		
		coeFourierLimited = new Grid1DComplex (momentCurveLimited, extendToNextPowerOfTwo);
		coeFourierLimited.transformForward();
		consistencyConstrain(coeFourierLimited, n);
		momentCurveRestored0 = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
		
		momentCurveRestored = momentCurveRestored0;
		float err = 0, val;
		for(int l = 0; l < maxIter; l++){

			for(int i = 0; i < numThetaLimited; i++){
				momentCurveRestored.setAtIndex(i, 0);
				momentCurveRestored.setAtIndex(i+numThetaPI, 0);
			}
			
			coeFourierLimited = new Grid1DComplex (momentCurveRestored, extendToNextPowerOfTwo);
			coeFourierLimited.transformForward();
			consistencyConstrain(coeFourierLimited, n);
			momentCurveRestored = getMomentCurveFromFourierCoefficients(coeFourierLimited, momentCurve.getNumberOfElements());
			err = 0;
			for(int i = 0; i < momentCurve.getNumberOfElements(); i++){
				momentCurveRestored.setAtIndex(i, momentCurveRestored.getAtIndex(i) * 0.5f + momentCurveRestored0.getAtIndex(i));
				val = momentCurveRestored.getAtIndex(i) - momentCurve.getAtIndex(i);
				err = err + val * val;
				}
			err = err/ (numThetaPI-numThetaLimited);
			err = (float)Math.sqrt(err);
			System.out.println("n = " + n + " error = " + err);
		}
	
		return momentCurveRestored;
	}
	/**
	 * get the matrix X when n is even, X is for complete sinogram
	 * @param n: the order which is even
	 * @param sinogram: complete sinogram 
	 * @return
	 */
	public SimpleMatrix getMatrixXEven(int n, Grid2D sinogram){
		int thetaNumTotal = sinogram.getSize()[1];
		double deltaTheta = sinogram.getSpacing()[1];
		double [][] matX = new double[thetaNumTotal][n + 1];
		for(int row = 0; row < thetaNumTotal; row ++)
		{
			matX[row][0] = 1;
			for(int col = 1; col <= n/2; col ++){
				matX[row][col * 2 - 1] = Math.cos(2 * col *deltaTheta * row);
				matX[row][col * 2 ] = Math.sin(2 * col * deltaTheta * row);
			}
		}
		SimpleMatrix X =new SimpleMatrix(new Matrix(matX));
		return X;
	}
	
	
	/**
	 * get matrix X when n is odd, X is for complete sinogram
	 * @param n: should be odd
	 * @param sinogram: complete sinogram
	 * @return
	 */
	public SimpleMatrix getMatrixXOdd(int n, Grid2D sinogram){
		int thetaNumTotal = sinogram.getSize()[1];
		double deltaTheta = sinogram.getSpacing()[1];
		double [][] matX = new double[thetaNumTotal][n + 1];
		for(int row = 0; row < thetaNumTotal; row ++)
		{
			for(int col = 0; col <= (n-1)/2; col ++){
				matX[row][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * row);
				matX[row][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * row);
			}
		}
		SimpleMatrix X =new SimpleMatrix(new Matrix(matX));
		return X;
	}
	
	/**
	 * get the matrix X for n from the matrix of n-2 to save computation, when n is even or odd
	 * @param n
	 * @param sinogram: complete sinogram
	 * @param Xpre: matrix X of order n - 2
	 * @return
	 */
	public SimpleMatrix getMatrixX(int n, Grid2D sinogram, SimpleMatrix Xpre){
		double deltaTheta = sinogram.getSpacing()[1];
		SimpleMatrix X = new SimpleMatrix(Xpre.getRows(),n+1);
		X.setSubMatrixValue(0, 0, Xpre);
		for(int row = 0; row < Xpre.getRows(); row ++){
			X.setElementValue(row, n - 1, Math.cos(n*deltaTheta*row));
			X.setElementValue(row, n, Math.sin(n*deltaTheta*row));
		}
		return X;
	}
	
	/**
	 * get the moment curves a_n(\theta) for the complete sinogram
	 * @param n
	 * @param sinogram
	 * @return
	 */
	public SimpleVector getVectorY(int n, Grid2D sinogram){
		Grid1D momentSum = ChebyshevFourierTransformMoment(sinogram, n);
		SimpleVector Y = new SimpleVector(momentSum.getBuffer());	
		return Y;
	}
	
/**
 * estimate the missing part of the moment curve
 * @param order
 * @param numThetaLimited: number of limited angle theta
 * @param Xcomplete: Matrix X when full data
 * @param Ycomplete: Vector Y when full data
 * @return
 */
	
	public SimpleVector getEstimatedMomentCurve(int order, int numThetaLimited, SimpleMatrix Xcomplete, SimpleVector Ycomplete){
		SimpleMatrix XLimited = Xcomplete.getSubMatrix(0, 0, numThetaLimited, Xcomplete.getCols());
		SimpleVector YLimited = Ycomplete.getSubVec(0, numThetaLimited);
		double tau = 0.0005 * (1 - order/1000); 
		iterativeSoftThresholding lasso = new iterativeSoftThresholding();
	
		SimpleVector beta =lasso.runLasso(XLimited, YLimited, tau);
		SimpleVector momentCurveRecovered = SimpleOperators.multiply(Xcomplete, beta);
		
		return momentCurveRecovered;
	}
	
	/**
	 * estimate the missing part of the moment curve using Weighted Least Square estimation (WLS)
	 * @param numThetaLimited
	 * @param Xcomplete
	 * @param Ycomplete
	 * @return
	 */
	public SimpleVector getEstimatedMomentCurveWLS( int order,int numThetaLimited, SimpleMatrix Xcomplete, SimpleVector Ycomplete){
		SimpleMatrix XLimited = Xcomplete.getSubMatrix(0, 0, numThetaLimited, Xcomplete.getCols());
		SimpleVector YLimited = Ycomplete.getSubVec(0, numThetaLimited);
	
		RobustFitting robust = new RobustFitting();
	
		SimpleVector beta =robust.iterativeReweightedLeastSquares(order,XLimited, YLimited);
		SimpleVector momentCurveRecovered = SimpleOperators.multiply(Xcomplete, beta);
		
		return momentCurveRecovered;
	}

	
	/**
	 * estimate the missing part of the moment curve using ridge regression
	 * @param numThetaLimited
	 * @param Xcomplete
	 * @param Ycomplete
	 * @return
	 */
	public SimpleVector getEstimatedMomentCurveRidge( float lambda,int numThetaLimited, SimpleMatrix Xcomplete, SimpleVector Ycomplete){
		SimpleMatrix XLimited = Xcomplete.getSubMatrix(0, 0, numThetaLimited, Xcomplete.getCols());
		SimpleVector YLimited = Ycomplete.getSubVec(0, numThetaLimited);
	
		RobustFitting robust = new RobustFitting();
	
		SimpleVector beta =robust.RidgeRegression(XLimited, YLimited, lambda);
		SimpleVector momentCurveRecovered = SimpleOperators.multiply(Xcomplete, beta);
		
		return momentCurveRecovered;
	}
	
	public Grid1D getSparsityOfMomentCurves(Grid2D momentCurveMatrixGT){
		Grid1D sparsity = new Grid1D(momentCurveMatrixGT.getHeight());
		for(int i = 0; i < sparsity.getNumberOfElements(); i++){
			sparsity.setAtIndex(i, this.getSparsityOfFourierParameters(i, momentCurveMatrixGT.getSubGrid(i)));
		}
		
		return sparsity;
	}
	
	/**
	 * get the number of parameters whose absolute values are bigger than the threshold
	 * @param order
	 * @param momentCurve
	 * @return
	 */
	public int getSparsityOfFourierParameters(int order,  Grid1D momentCurve){
		
		//Grid1D momentCurve = this.vector2Grid(curve);
		Grid1DComplex coe = this.ChebyshevFourierTransformFromMomentCurve(momentCurve);
		double tau = 0.001 * (1 - order/1000); //threshold in Lasso regression
		double tau2 = tau * momentCurve.getNumberOfElements()/2; //relation between Fourier coefficients and the cosine and sine basis parameters
		int count = 0;
		if(order%2==0){//even cases	
			if(Math.abs(coe.getRealAtIndex(0)) > (tau2*2)) //DC component
				count++;
			for(int i = 2; i <= order; i=i+2)
			{
				if(Math.abs(coe.getRealAtIndex(i)) > tau2) 
					count++;
				if(Math.abs(coe.getImagAtIndex(i)) > tau2)
					count++;
			}
		}
		else//odd cases
		{
			for(int i = 1; i <= order; i=i+2)
			{
				if(Math.abs(coe.getRealAtIndex(i)) > tau2) 
					count++;
				if(Math.abs(coe.getImagAtIndex(i)) > tau2)
					count++;
			}
		}
		

		return count;
	}
	
	/**
	 * 
	 * @param vec
	 * @return
	 */
	
	private Grid1D vector2Grid(SimpleVector vec){
		Grid1D grid = new Grid1D(vec.getLen());
		for(int i = 0; i < vec.getLen(); i++)
			grid.setAtIndex(i, (float) (vec.getElement(i)));
		
		return grid;
	}

	/**
	 * Lasso regression to estimate the parameters
	 * @param X
	 * @param Y
	 * @param n
	 * @return
	 */
	public SimpleVector ChebyshevFourierTransformCoefficientLassoRegression(SimpleMatrix X, SimpleVector Y, double tau){
		//Normalize to guarantee convergence		
		iterativeSoftThresholding lasso = new iterativeSoftThresholding();
		//double tau = 0.001 * (1 - n/1000); // set tau
		SimpleVector beta = lasso.runLasso(X, Y, tau);
		return beta;
	}
	
	/**
	 * calculate the nth moment of the limited angle sinogram, i.e. angles 0 - 160, based on the checkboard pattern of 
	 * the Chebyshev Fourier transform of the ideal complete sinogram, we can get the moment values for the angles 161 - 180 
	 * using Fourier curve fitting algorithm, which can be simplified to be a multiple linear regression problem
	 * @param sinogram  limited angle sinogram
	 * @param n   nth order moment
	 * @return
	 */
	public Grid1DComplex ChebyshevFourierTransformCoefficientRegression(Grid2D sinogram, int n){
		Grid1D momentSum = ChebyshevFourierTransformMoment(sinogram, n);
		
		int thetaNum = sinogram.getSize()[1];
		int thetaNumTotal = (int) (2 * Math.PI / sinogram.getSpacing()[1]);// the number of projections in 360 full scan
		
		Grid1DComplex coeFourier = new Grid1DComplex (new Grid1D(thetaNumTotal), false);
		double deltaTheta = sinogram.getSpacing()[1];
		
		if(n < thetaNum){
		if(n == 0){
			float val = (float) (momentSum.getGridOperator().sum(momentSum)/thetaNum * thetaNumTotal);
			coeFourier.setRealAtIndex(0, val);
		}
		else {
		    
			double [][] matX = new double[thetaNum][n + 1];
			double [][] vecY = new double[thetaNum][1];
			if (n % 2 == 0){
			/*
			 * [1 cos(2*deltaTheta) sin(2*deltaTheta) cos(4*deltaTheta) cos(4 * deltaTheta) ... cos(n*deltaTheta) sin(n*deltaTheta)
			 * ...
			 * 1 cos(2*deltaTheta*rowIdx) sin(2*deltaTheta*rowIdx) cos(4*deltaTheta*rowIdx) cos(4 * deltaTheta*rowIdx) ... cos(n*deltaTheta*rowIdx) sin(n*deltaTheta*rowIdx)
			 * ...]
			 */
			
			for(int row = 0; row < thetaNum; row ++)
			{
				vecY[row][0] = momentSum.getAtIndex(row);
				matX[row][0] = 1;
				for(int col = 1; col <= n/2; col ++){
					matX[row][col * 2 - 1] = Math.cos(2 * col *deltaTheta * row);
					matX[row][col * 2 ] = Math.sin(2 * col * deltaTheta * row);
				}
			}
			
		}
		else{
			/*
			 * [ cos(deltaTheta) sin(deltaTheta) cos(3*deltaTheta) cos(3 * deltaTheta) ... cos(n*deltaTheta) sin(n*deltaTheta)
			 * ...
			 *  cos(deltaTheta*rowIdx) sin(deltaTheta*rowIdx) cos(3*deltaTheta*rowIdx) cos(3 * deltaTheta*rowIdx) ... cos(n*deltaTheta*rowIdx) sin(n*deltaTheta*rowIdx)
			 * ...]
			 */		
			for(int row = 0; row < thetaNum; row ++)
			{
				vecY[row][0] = momentSum.getAtIndex(row);
				for(int col = 0; col <= (n-1)/2; col ++){
					matX[row][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * row);
					matX[row][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * row);
				}
			}				
		}
			SimpleMatrix X =new SimpleMatrix(new Matrix(matX));
			SimpleVector Y = new SimpleVector(new Matrix(vecY));
			SimpleMatrix X_pinv = X.inverse(SimpleMatrix.InversionType.INVERT_SVD);
			SimpleVector vecCoe = SimpleOperators.multiply(X_pinv, Y);
			
			if (n % 2 == 0){
				coeFourier.setRealAtIndex(0, (float)(vecCoe.getElement(0) * thetaNumTotal));
				for(int idx = 1; idx <= n/2; idx ++){
				coeFourier.setRealAtIndex(idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(idx * 2, (float) (- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
				coeFourier.setRealAtIndex(thetaNumTotal - idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(thetaNumTotal - idx * 2, (float)(- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));		
				}
			}else{
				for(int idx = 0; idx <= (n-1)/2; idx ++){
					coeFourier.setRealAtIndex(idx * 2 + 1, (float)(vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(idx * 2 + 1, (float) (- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));
					coeFourier.setRealAtIndex(thetaNumTotal - idx * 2 - 1, (float)(vecCoe.getElement(idx * 2) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(thetaNumTotal - idx * 2 - 1, (float)(- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));		
					}
			}
		}
		}
		else
			System.out.println("the order n is too large! n should be smaller than thetaNum!");
								
		return coeFourier;
	}
	
	
	/**
	 * calculate the nth moment of the limited angle sinogram, i.e. angles 0 - 160, based on the checkboard pattern of 
	 * the Chebyshev Fourier transform of the ideal complete sinogram, we can get the moment values for the angles 161 - 180 
	 * using Fourier curve fitting algorithm, which can be simplified to be a multiple linear regression problem
	 * @param sinogram  limited angle sinogram, but with full size, i.e. 180 - 340 degree data are added based on symmetric property
	 * @param n   nth order moment
	 * @param thetaNum the number of projections in limited angle sinogram
	 * @return
	 */
	public Grid1DComplex ChebyshevFourierTransformCoefficientRegression2(Grid2D sinogram, int n, int thetaNum){
		Grid1D momentSum = ChebyshevFourierTransformMoment(sinogram, n);
		int thetaNumTotal = (int) (2 * Math.PI / sinogram.getSpacing()[1]);// the number of projections in 360 full scan
		int halfNum = (int) (thetaNumTotal/2);
		Grid1DComplex coeFourier = new Grid1DComplex (new Grid1D(thetaNumTotal), false);
		double deltaTheta = sinogram.getSpacing()[1];
		
		if(n < thetaNum * 2){
		if(n == 0){
			float val = 0;
			for (int thetaIdx = 0; thetaIdx < thetaNum; thetaIdx++){
			val = val + momentSum.getAtIndex(thetaIdx) + momentSum.getAtIndex(thetaIdx + halfNum);
			}
			val = val/(2 * thetaNum) * thetaNumTotal;
			coeFourier.setRealAtIndex(0, val);
		}
		else {
		    
			double [][] matX = new double[thetaNum * 2][n + 1];
			double [][] vecY = new double[thetaNum * 2][1];
			if (n % 2 == 0){
			/*
			 * [1 cos(2*deltaTheta) sin(2*deltaTheta) cos(4*deltaTheta) cos(4 * deltaTheta) ... cos(n*deltaTheta) sin(n*deltaTheta)
			 * ...
			 * 1 cos(2*deltaTheta*rowIdx) sin(2*deltaTheta*rowIdx) cos(4*deltaTheta*rowIdx) cos(4 * deltaTheta*rowIdx) ... cos(n*deltaTheta*rowIdx) sin(n*deltaTheta*rowIdx)
			 * ...]
			 */
			
			for(int row = 0; row < thetaNum; row ++)
			{
				vecY[row][0] = momentSum.getAtIndex(row);
				vecY[row + thetaNum][0] = momentSum.getAtIndex(row + halfNum);
				matX[row][0] = 1;
				matX[row + thetaNum][0] = 1;
				for(int col = 1; col <= n/2; col ++){
					matX[row][col * 2 - 1] = Math.cos(2 * col *deltaTheta * row);
					matX[row][col * 2 ] = Math.sin(2 * col * deltaTheta * row);
					matX[row + thetaNum][col * 2 - 1] = Math.cos(2 * col *deltaTheta * (row + halfNum));
					matX[row + thetaNum][col * 2 ] = Math.sin(2 * col * deltaTheta * (row + halfNum));
				}
			}
			
		}
		else{
			/*
			 * [ cos(deltaTheta) sin(deltaTheta) cos(3*deltaTheta) cos(3 * deltaTheta) ... cos(n*deltaTheta) sin(n*deltaTheta)
			 * ...
			 *  cos(deltaTheta*rowIdx) sin(deltaTheta*rowIdx) cos(3*deltaTheta*rowIdx) cos(3 * deltaTheta*rowIdx) ... cos(n*deltaTheta*rowIdx) sin(n*deltaTheta*rowIdx)
			 * ...]
			 */		
			for(int row = 0; row < thetaNum; row ++)
			{
				vecY[row][0] = momentSum.getAtIndex(row);
				vecY[row + thetaNum][0] = momentSum.getAtIndex(row + halfNum);
				for(int col = 0; col <= (n-1)/2; col ++){
					matX[row][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * row);
					matX[row][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * row);
					matX[row + thetaNum][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * (row + halfNum));
					matX[row + thetaNum][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * (row + halfNum));
				}
			}				
		}
			SimpleMatrix X =new SimpleMatrix(new Matrix(matX));
			SimpleVector Y = new SimpleVector(new Matrix(vecY));
			SimpleMatrix X_pinv = X.inverse(SimpleMatrix.InversionType.INVERT_SVD);
			SimpleVector vecCoe = SimpleOperators.multiply(X_pinv, Y);
			
			if (n % 2 == 0){
				coeFourier.setRealAtIndex(0, (float)(vecCoe.getElement(0) * thetaNumTotal));
				for(int idx = 1; idx <= n/2; idx ++){
				coeFourier.setRealAtIndex(idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(idx * 2, (float) (- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
				coeFourier.setRealAtIndex(thetaNumTotal - idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(thetaNumTotal - idx * 2, (float)(- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));		
				}
			}else{
				for(int idx = 0; idx <= (n-1)/2; idx ++){
					coeFourier.setRealAtIndex(idx * 2 + 1, (float)(vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(idx * 2 + 1, (float) (- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));
					coeFourier.setRealAtIndex(thetaNumTotal - idx * 2 - 1, (float)(vecCoe.getElement(idx * 2) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(thetaNumTotal - idx * 2 - 1, (float)(- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));		
					}
			}
		}
		}
		else
			System.out.println("the order n is too large! n should be smaller than thetaNum!");
								
		return coeFourier;
	}
	
	public Grid1DComplex ChebyshevFourierTransformCoefficientRegression3(Grid2D sinogram, int n, int thetaNum, int m){
		Grid1D momentSum = ChebyshevFourierTransformMoment(sinogram, n);
		int thetaNumTotal = (int) (2 * Math.PI / sinogram.getSpacing()[1]);// the number of projections in 360 full scan
		int halfNum = (int) (thetaNumTotal/2);
		Grid1DComplex coeFourier = new Grid1DComplex (new Grid1D(thetaNumTotal), false);
		double deltaTheta = sinogram.getSpacing()[1];
		
		if(n < thetaNum * 2){//not under-determined
		if(n == 0){
			float val = 0;
			for (int thetaIdx = 0; thetaIdx < thetaNum; thetaIdx++){
			val = val + momentSum.getAtIndex(thetaIdx) + momentSum.getAtIndex(thetaIdx + halfNum);
			}
			val = val/(2 * thetaNum) * thetaNumTotal;
			coeFourier.setRealAtIndex(0, val);
		}
		else {
		    
			double [][] matX = new double[thetaNum * 2][m + 1];
			double [][] vecY = new double[thetaNum * 2][1];
			if (n % 2 == 0){
			/*
			 * [1 cos(2*deltaTheta) sin(2*deltaTheta) cos(4*deltaTheta) cos(4 * deltaTheta) ... cos(m*deltaTheta) sin(m*deltaTheta)
			 * ...
			 * 1 cos(2*deltaTheta*rowIdx) sin(2*deltaTheta*rowIdx) cos(4*deltaTheta*rowIdx) cos(4 * deltaTheta*rowIdx) ... cos(m*deltaTheta*rowIdx) sin(m*deltaTheta*rowIdx)
			 * ...]
			 */
			
			for(int row = 0; row < thetaNum; row ++)
			{
				vecY[row][0] = momentSum.getAtIndex(row);
				vecY[row + thetaNum][0] = momentSum.getAtIndex(row + halfNum);
				matX[row][0] = 1;
				matX[row + thetaNum][0] = 1;
				for(int col = 1; col <= m/2; col ++){
					matX[row][col * 2 - 1] = Math.cos(2 * col *deltaTheta * row);
					matX[row][col * 2 ] = Math.sin(2 * col * deltaTheta * row);
					matX[row + thetaNum][col * 2 - 1] = Math.cos(2 * col *deltaTheta * (row + halfNum));
					matX[row + thetaNum][col * 2 ] = Math.sin(2 * col * deltaTheta * (row + halfNum));
				}
			}
			
			}
			else{
			/*
			 * [ cos(deltaTheta) sin(deltaTheta) cos(3*deltaTheta) cos(3 * deltaTheta) ... cos(m*deltaTheta) sin(m*deltaTheta)
			 * ...
			 *  cos(deltaTheta*rowIdx) sin(deltaTheta*rowIdx) cos(3*deltaTheta*rowIdx) cos(3 * deltaTheta*rowIdx) ... cos(m*deltaTheta*rowIdx) sin(m*deltaTheta*rowIdx)
			 * ...]
			 */		
				for(int row = 0; row < thetaNum; row ++)
				{
					vecY[row][0] = momentSum.getAtIndex(row);
					vecY[row + thetaNum][0] = momentSum.getAtIndex(row + halfNum);
					for(int col = 0; col <= (m-1)/2; col ++){
						matX[row][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * row);
						matX[row][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * row);
						matX[row + thetaNum][col * 2 ] = Math.cos((2 * col + 1) *deltaTheta * (row + halfNum));
						matX[row + thetaNum][col * 2 + 1] = Math.sin((2 * col + 1) * deltaTheta * (row + halfNum));
					}
				}				
			}
			SimpleMatrix X =new SimpleMatrix(new Matrix(matX));
			SimpleVector Y = new SimpleVector(new Matrix(vecY));
			
			RobustFitting rft = new RobustFitting();
			//SimpleMatrix X_pinv = X.inverse(SimpleMatrix.InversionType.INVERT_SVD);
			//SimpleVector vecCoe = SimpleOperators.multiply(X_pinv, Y);
			//SimpleVector vecCoe =rft.ordinaryLeastSquares(X, Y);
			//SimpleVector vecCoe =rft.ordinaryLeastSquaresPseudoInverse(X, Y);
			SimpleVector vecCoe = rft.weightedLeastSquares(X, Y, 0.1);
			
			if (n % 2 == 0){
				coeFourier.setRealAtIndex(0, (float)(vecCoe.getElement(0) * thetaNumTotal));
				for(int idx = 1; idx <= m/2; idx ++){
				coeFourier.setRealAtIndex(idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(idx * 2, (float) (- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
				coeFourier.setRealAtIndex(thetaNumTotal - idx * 2, (float)(vecCoe.getElement(idx * 2 -1) * thetaNumTotal / 2));
				coeFourier.setImagAtIndex(thetaNumTotal - idx * 2, (float)(- vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));		
				}
			}else{
				for(int idx = 0; idx <= (m-1)/2; idx ++){
					coeFourier.setRealAtIndex(idx * 2 + 1, (float)(vecCoe.getElement(idx * 2 ) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(idx * 2 + 1, (float) (- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));
					coeFourier.setRealAtIndex(thetaNumTotal - idx * 2 - 1, (float)(vecCoe.getElement(idx * 2) * thetaNumTotal / 2));
					coeFourier.setImagAtIndex(thetaNumTotal - idx * 2 - 1, (float)(- vecCoe.getElement(idx * 2 + 1) * thetaNumTotal / 2));		
					}
			}
		}
		}
		else
			System.out.println("the order n is too large! n should be smaller than thetaNum!");
								
		return coeFourier;
	}
	
	/**
	 * when the object support is between [-r, r] instead of [-1, 1], where 0<r<1, the triangle region becomes thinner
	 * the slope is about r instead of 1
	 * @param maxOrder
	 * @param numDetector
	 * @param r
	 * @param thres
	 * @return
	 */
	public Grid1D getHelgasonLudwigSlope(int maxOrder, int numDetector, float r, float thres ){
		double deltaS = 1./numDetector;
		double s = -1 + 0.5 *deltaS;
		Grid1D sPowerk = new Grid1D(numDetector);
		Grid1D ss = new Grid1D(numDetector);
		Grid1D rss = new Grid1D(numDetector);
		Grid1D tss = new Grid1D(numDetector);
		
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			s = -1 + (i + 0.5) * deltaS;
			//s = Math.sin((n + 1) * Math.acos(s))/Math.sin(Math.acos(s)); //second kind Chebyshev
			ss.setAtIndex(i, (float) s); 
			rss.setAtIndex(i, (float) (s * r));
			tss.setAtIndex(i, (float) (Math.sqrt(1 - s * s)));
		}
		
		Grid1D slopeN = new Grid1D(maxOrder + 1);
		
		Grid1D sPowerk2 = new Grid1D(numDetector);
		double val;
		boolean tag = false;
		for(int n2 = 0; n2 <= maxOrder; n2++){
			for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
				s = rss.getAtIndex(i);
				s = Math.sin((n2 + 1) * Math.acos(s))/Math.sin(Math.acos(s)); //second kind Chebyshev
				sPowerk2.setAtIndex(i, (float) s);  
			}
			tag = false;
			
			for(int m = n2; m >= 0; m = m - 2){
				for(int n = m; n <= n2; n = n + 2){
					for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
						s = ss.getAtIndex(i);
						s = Math.sin((n + 1) * Math.acos(s))/Math.sin(Math.acos(s)); //second kind Chebyshev
						sPowerk.setAtIndex(i, (float) (s * tss.getAtIndex(i) * sPowerk2.getAtIndex(i))); 
					}
					val = sPowerk.getGridOperator().sum(sPowerk);
					if(Math.abs(val) > thres){
						tag = true;
						break;
					}
						
				}
				if(tag)
				{
					//if(debug)
					//System.out.println(n2 + " " + m);
					slopeN.setAtIndex(n2, m);
					break;
				}
			}
			
		}
		
		return slopeN;
	}
	
	public Grid1D getSlope(int maxOrder){
		double[] slope1 = new double[]{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 26.0, 27.0, 28.0, 29.0, 30.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 41.0, 42.0, 43.0, 44.0, 45.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 48.0, 49.0, 50.0, 51.0, 52.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 55.0, 56.0, 57.0, 58.0, 59.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 62.0, 63.0, 64.0, 65.0, 66.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 69.0, 70.0, 71.0, 72.0, 73.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 76.0, 77.0, 78.0, 79.0, 80.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 83.0, 84.0, 85.0, 86.0, 87.0, 86.0, 87.0, 88.0, 89.0, 90.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 93.0, 94.0, 95.0, 96.0, 97.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 100.0, 101.0, 102.0, 103.0, 104.0, 103.0, 104.0, 105.0, 106.0, 107.0, 106.0, 107.0, 108.0, 109.0, 110.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 114.0, 115.0, 116.0, 117.0, 118.0, 117.0, 118.0, 119.0, 120.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 127.0, 128.0, 129.0, 130.0, 131.0, 130.0, 131.0, 132.0, 133.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 138.0, 139.0, 140.0, 141.0, 142.0, 141.0, 142.0, 143.0, 144.0, 143.0, 144.0, 145.0, 146.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 151.0, 152.0, 153.0, 154.0, 155.0, 154.0, 155.0, 156.0, 157.0, 156.0, 157.0, 158.0, 159.0, 160.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 164.0, 165.0, 166.0, 167.0, 168.0, 167.0, 168.0, 169.0, 170.0, 169.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 174.0, 175.0, 176.0, 177.0, 178.0, 177.0, 178.0, 179.0, 180.0, 181.0, 180.0, 181.0, 182.0, 183.0, 184.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 187.0, 188.0, 189.0, 190.0, 189.0, 190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 194.0, 195.0, 196.0, 197.0, 198.0, 197.0, 198.0, 199.0, 200.0, 199.0, 200.0, 201.0, 202.0, 203.0, 202.0, 203.0, 206.0, 205.0, 206.0, 207.0, 208.0, 209.0, 208.0, 209.0, 210.0, 211.0, 210.0, 211.0, 212.0, 213.0, 212.0, 213.0, 214.0, 215.0, 214.0, 215.0, 218.0, 219.0, 220.0, 219.0, 220.0, 221.0, 222.0, 221.0, 222.0, 223.0, 222.0, 223.0, 224.0, 225.0, 224.0, 225.0, 226.0, 227.0, 230.0, 231.0, 230.0, 231.0, 232.0};
		Grid1D slope = new Grid1D(maxOrder);
		for (int n = 0; n <maxOrder; n++){
			if(n<361)
				slope.setAtIndex(n, (float)(slope1[n]));
			else
				slope.setAtIndex(n, (float)(Math.floor(n*0.6146)));
		}
		
		return slope;
	}
	
	
	/**
	 * get nth order ChebyshevFourierTransform Moment curves a_n(\theta), before Fourier transform
	 * @param sinogram
	 * @param n
	 * @return
	 */
public Grid1D ChebyshevFourierTransformMoment(Grid2D sinogram, int n){
		
		Grid1D momentSum = new Grid1D(sinogram.getSize()[1]);
		double s;
		double deltaS = 2.f / sinogram.getSize()[0]; // here regard the detector from -1 to 1, maxS = 2, deltaS = 2 / detectorSize;
		
		Grid1D sPowerk = new Grid1D(sinogram.getSize()[0]);
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			s = -1 + (i + 0.5) * deltaS;
			//s = Math.cos(n * Math.acos(s)); //Chebyshev first kind
			s = Math.sin((n + 1) * Math.acos(s))/Math.sin(Math.acos(s)); //second kind Chebyshev
			sPowerk.setAtIndex(i, (float) s);  //FIXME
		}
			
		Grid1D sub;
		for ( int i = 0; i < sinogram.getSize()[1]; i ++){
			sub = new Grid1D(sinogram.getSubGrid(i));
			sPowerk.getGridOperator().multiplyBy(sub, sPowerk);
			momentSum.setAtIndex(i, (float) (deltaS * sub.getGridOperator().sum(sub)));
		}
	
	
		return momentSum;
	}
	
/**
 * return the magnitudes of the ChebyshevFourier transform coefficients
 * @param sinogram
 * @param n
 * @return
 */
	public Grid1D ChebyshevFourierTransformMag(Grid2D sinogram, int n){
		Grid1D coefs;
		coefs= ChebyshevFourierTransform(sinogram, n).getMagSubGrid(0, sinogram.getHeight());
		return coefs;
	}
	
	
	/**
	 * add nth order moment components based on the Chebyshev-Fourier data consistency
	 * @param coeFourier
	 * @param n
	 * @param sinogram
	 */
	public void inverseChebyshevFourierTransform(Grid1DComplex coeFourier, int n, Grid2D sinogram){
		
	
		coeFourier.transformInverse();
		
		Grid1D coes = coeFourier.getRealSubGrid(0, sinogram.getSize()[1]);
	
		inverseChebyshevFourierTransformFromMomentCurves(coes, n, sinogram);
		
	}
	
	public Grid1D getMomentCurveFromFourierCoefficients(Grid1DComplex coeFourier, int length){
	
		coeFourier.transformInverse();
		
		Grid1D curve = coeFourier.getRealSubGrid(0, length);
		return curve;
	}
	
	public void inverseChebyshevFourierTransformFromMomentCurves(Grid1D curve, int n, Grid2D sinogram){
		double s, Ts, as;
		
		double deltaS = 2.f / sinogram.getSize()[0]; // here regard the detector from -1 to 1, maxS = 2, deltaS = 2 / detectorSize;
		
		Grid1D sPowerk = new Grid1D(sinogram.getSize()[0]);
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			s = -1 + (i + 0.5) * deltaS;
			as = Math.sqrt(1 - s * s);
			//Ts = Math.cos(n * Math.acos(s));
			Ts = Math.sin((n + 1) * Math.acos(s))/Math.sin(Math.acos(s));
			sPowerk.setAtIndex(i, (float) (as * Ts * 2.0/ Math.PI));  //FIXME
		}
		
		float val;
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			for(int j = 0; j < sinogram.getSize()[1]; j ++){
				val = (float)( sinogram.getAtIndex(i, j) + curve.getAtIndex(j) * sPowerk.getAtIndex(i));
				sinogram.setAtIndex(i, j, val);
			}
		}
		
	}
	
	
	/**
	 * add nth order moment components based on the Chebyshev-Fourier data consistency
	 * @param coeFourier
	 * @param n
	 * @param coe2DComplex  
	 */
	public void inverseChebyshevFourierTransformStep1(Grid1DComplex coeFourier, int n, Grid2DComplex coe2DComplex){		
		double s, Ts, as;

		double deltaS = 2.f / coe2DComplex.getSize()[0]; // here regard the detector from -1 to 1, maxS = 2, deltaS = 2 / detectorSize;
		
		Grid1D sPowerk = new Grid1D(coe2DComplex.getSize()[0]);

		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			s = -1 + (i + 0.5) * deltaS;
			as = Math.sqrt(1 - s * s);
			//Ts = Math.cos(n * Math.acos(s));
			Ts = Math.sin((n + 1) * Math.acos(s))/Math.sin(Math.acos(s));
			sPowerk.setAtIndex(i, (float) (as * Ts * 2.0/ Math.PI));  //FIXME
		}
		
	
		for( int i = 0; i < sPowerk.getNumberOfElements(); i++){
			for(int j = 0; j < coe2DComplex.getSize()[1]; j ++){
			coe2DComplex.setRealAtIndex(i, j, coe2DComplex.getRealAtIndex(i, j) + coeFourier.getRealAtIndex(j) * sPowerk.getAtIndex(i));
			coe2DComplex.setImagAtIndex(i, j, coe2DComplex.getImagAtIndex(i, j) + coeFourier.getImagAtIndex(j) * sPowerk.getAtIndex(i));
			}
		}
	}
	
	
	public void inverseChebyshevFourierTransformStep2(Grid2DComplex coe2DComplex, Grid2D sinogram){
		Grid1DComplex coe1DComplex = new Grid1DComplex(new Grid1D(sinogram.getSize()[1]));
		Grid1D sinoColumn;
		for(int i = 0; i < sinogram.getSize()[0]; i ++){
			for(int j = 0; j < sinogram.getSize()[1]; j ++){
				coe1DComplex.setRealAtIndex(j, coe2DComplex.getRealAtIndex(i, j));
				coe1DComplex.setImagAtIndex(j, coe2DComplex.getImagAtIndex(i, j));
			}
			coe1DComplex.transformInverse();
			sinoColumn = coe1DComplex.getRealSubGrid(0, sinogram.getSize()[1]);
			for(int j = 0; j < sinogram.getSize()[1]; j ++){
				sinogram.setAtIndex(i, j, sinoColumn.getAtIndex(j));
			}
		}
	}

	
	/**
	 * satisfy the Helgason-Ludwig data consistency , the check board pattern
	 * here coeFourier is not fftshifted
	 * @param coeFourier
	 * @param n
	 */
	public void consistencyConstrain(Grid1DComplex coeFourier, int n) {

		int L = coeFourier.getSize()[0];
		if (2 * n + 1 <= L) {
			if (n % 2 == 0) {
				for (int i = 1; i <= n; i = i + 2) {
					coeFourier.setImagAtIndex(i, 0);
					coeFourier.setRealAtIndex(i, 0);
				}
				for (int i = n + 1; i < L - n; i++) {
					coeFourier.setImagAtIndex(i, 0);
					coeFourier.setRealAtIndex(i, 0);
				}
				if (L % 2 == 0)
					for (int i = L - n + 1; i < L; i = i + 2) {
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
				else
					for (int i = L - n + 2; i < L; i = i + 2) {
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
			} else {
				for (int i = 0; i <= n; i = i + 2) {
					coeFourier.setImagAtIndex(i, 0);
					coeFourier.setRealAtIndex(i, 0);
				}
				for (int i = n + 1; i < L - n; i++) {
					coeFourier.setImagAtIndex(i, 0);
					coeFourier.setRealAtIndex(i, 0);
				}
				if (L % 2 == 0)
					for (int i = L - n + 1; i < L; i = i + 2) {
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
				else
					for (int i = L - n + 2; i < L; i = i + 2) {
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
			}
		}
		else{
			if(L % 2 == 0){
				if(n % 2 == 0)
					for(int i = 1; i < L; i = i + 2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
				else
					for(int i = 0; i < L; i = i + 2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
			}
			else{
				if(n % 2 == 0){
					for(int i = 1; i <= (L-1)/2; i = i + 2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
					for(int i = L -1; i >= (L+1)/2; i = i -2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
				}
				else{
					for(int i = 0; i <= (L-1)/2; i = i + 2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
					for(int i = L -2; i >= (L+1)/2; i = i -2){
						coeFourier.setImagAtIndex(i, 0);
						coeFourier.setRealAtIndex(i, 0);
					}
				}
			}
						
		}
	}
	
	public void fftshift1D(Grid1D magnFourier){
		Grid1D copy = new Grid1D(magnFourier);
		int len = magnFourier.getNumberOfElements();
		int halfLen;
		if(len % 2 == 0)
			 halfLen = (int) Math.ceil(len/2);
		else
			halfLen = (int) Math.ceil(len/2 + 1);
		
		for (int i = 0; i < len - halfLen; i ++)
			magnFourier.setAtIndex(i, copy.getAtIndex(i + halfLen));
		for(int i = len - halfLen; i < len; i ++)
			magnFourier.setAtIndex(i, copy.getAtIndex( i - (len - halfLen)));			
	}
	
	
	
	/**
	 *  cross-correlation coefficient r of two vectors
	 * @param curve1
	 * @param curve2
	 * @return
	 */
	public double correlationCoefficient(Grid1D curve1, Grid1D curve2){
		double mean1 = curve1.getGridOperator().sum(curve1)/curve1.getNumberOfElements();
		double mean2 = curve2.getGridOperator().sum(curve2)/curve2.getNumberOfElements();
		curve1.getGridOperator().subtractBy(curve1, (float)(mean1));
		curve2.getGridOperator().subtractBy(curve2, (float)(mean2));
		Grid1D cross = (Grid1D) curve1.clone();
		cross.getGridOperator().multiplyBy(cross, curve2);
		curve1.getGridOperator().multiplyBy(curve1, curve1); //square
		curve2.getGridOperator().multiplyBy(curve2, curve2);  
		double a = Math.sqrt(curve1.getGridOperator().sum(curve1)*curve2.getGridOperator().sum(curve2));
		double r;
		if(a == 0)
			r = 1;
		else
			r = cross.getGridOperator().sum(cross)/a;
		return r;
	}
	
	public Grid1D getCorrCoeBetweenEstimationAndGroundTruth(Grid2D curvesGT, Grid2D curvesRecovered, ParallelSinogramAndReconstruction  objSino){
		Grid1D correlation = new Grid1D(curvesGT.getSize()[1]);
		int numThetaLimited = (int)(objSino.maxThetaLimited/objSino.deltaTheta);
		int numThetaPI = (int)(Math.PI/objSino.deltaTheta);
		int lengthMiss = numThetaPI - numThetaLimited -1;
		double r;
		Grid1D singleCurveRecovered,singleCurveGT;
	for(int order = 0; order < curvesGT.getSize()[1]; order++){
		singleCurveRecovered = (Grid1D)curvesRecovered.getSubGrid(order).clone();
		singleCurveGT = (Grid1D) curvesGT.getSubGrid(order).clone();
		r = this.correlationCoefficient(singleCurveRecovered.getSubGrid(numThetaLimited+1, lengthMiss), singleCurveGT.getSubGrid(numThetaLimited+1, lengthMiss));
		correlation.setAtIndex(order,(float) (r));
		
	}	
		return correlation;
	}
}
