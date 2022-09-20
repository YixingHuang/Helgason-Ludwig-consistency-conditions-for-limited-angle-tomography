package edu.stanford.rsl.tutorial.hlcc;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.numerics.DecompositionQR;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class RobustFitting {
	SimpleMatrix matX;
	SimpleVector vecY;
	public int maxIter;
	public float eps;
	public double lambda = 0.1;
	public boolean isRegular = true;
	public boolean debug = false;
	
	public RobustFitting(){
		
	}
	
	public static enum WeightType {
		andrews,
		bisquare,
		cauchy,
		fair,
		huber,
		logistic		
	}
	
	public SimpleVector ordinaryLeastSquaresPseudoInverse(SimpleMatrix matX, SimpleVector vecY){
		SimpleMatrix X_pinv = matX.inverse(SimpleMatrix.InversionType.INVERT_SVD);
		SimpleVector vecCoe = SimpleOperators.multiply(X_pinv, vecY);
		return vecCoe;
	}
	
	public SimpleVector ordinaryLeastSquares(SimpleMatrix matX, SimpleVector vecY){
		SimpleMatrix matXT = matX.transposed();
		SimpleMatrix matXTX = SimpleOperators.multiplyMatrixProd(matXT, matX);
		SimpleMatrix matXTXinv = matXTX.inverse(SimpleMatrix.InversionType.INVERT_QR);
		SimpleMatrix Xpinv = SimpleOperators.multiplyMatrixProd(matXTXinv, matXT);
		SimpleVector vecCoe = SimpleOperators.multiply(Xpinv, vecY);
		return vecCoe;
	}
	
	
	public SimpleVector RidgeRegression(SimpleMatrix matX, SimpleVector vecY, double lambda){
		SimpleMatrix matXT = matX.transposed();
		SimpleMatrix matXTX = SimpleOperators.multiplyMatrixProd(matXT, matX);
		SimpleMatrix eye = new SimpleMatrix(matX.getCols(), matX.getCols());
		eye.identity();
		eye.multiplyBy(lambda);
		SimpleMatrix matXTXLambda = SimpleOperators.add(eye, matXTX);
		SimpleMatrix matXTXinv = matXTXLambda.inverse(SimpleMatrix.InversionType.INVERT_QR);
		SimpleMatrix Xpinv = SimpleOperators.multiplyMatrixProd(matXTXinv, matXT);
		SimpleVector vecCoe = SimpleOperators.multiply(Xpinv, vecY);
		return vecCoe;
	}
	
	public SimpleVector LassoRegression(SimpleMatrix matX, SimpleVector vecY, double tau){
		iterativeSoftThresholding lasso = new iterativeSoftThresholding();
		SimpleVector beta = lasso.runLasso(matX, vecY, tau);
		return beta;
	}
	
	
	public SimpleVector weightedLeastSquares(SimpleMatrix matX, SimpleVector vecY, double tune){
		SimpleVector vecYhat;
		SimpleVector resd, diff;
		double sigma;
		double w, rmse = 1;
		SimpleMatrix matWX = new SimpleMatrix(matX);
		SimpleVector vecWY = new SimpleVector(vecY);
		SimpleVector vecCoeNew;
		int maxIter = 50, iter = 0;
		
		SimpleVector vecCoe = ordinaryLeastSquares(matX, vecY); //first estimate
		if(debug)
	    	System.out.println(matX.getRows() + " " + matX.getCols());
		SimpleVector hfac = leverageFactor(matX);
		//System.out.println(hfac.toString());
		
	  while(iter < maxIter && rmse > 1.0e-5){
		vecYhat = SimpleOperators.multiply(matX, vecCoe);
	    resd = SimpleOperators.subtract(vecY, vecYhat);
	    //sigma = getSigma(resd);
	    sigma = madSigma(resd, matX.getCols());
	    if(debug)
	    	System.out.println("sigma = " + sigma);
	    sigma = Math.max(sigma, 1 );
		resd.divideBy(tune * sigma);
		resd.multiplyElementWiseBy(hfac);
		
		for(int i = 0; i < vecY.getLen(); i ++){
			w = fair(resd.getElement(i));
			//w = huber(resd.getElement(i));
			vecWY.setElementValue(i, w * vecY.getElement(i));
			for(int j = 0; j < matX.getCols(); j ++)
				matWX.setElementValue(i, j, w * matX.getElement(i, j));
		}
		
		if(this.isRegular)
			vecCoeNew = RidgeRegression(matWX, vecWY, this.lambda);
		else
			vecCoeNew = ordinaryLeastSquares(matWX, vecWY);
		//vecCoeNew = ordinaryLeastSquaresPseudoInverse(matWX, vecWY);
		diff = SimpleOperators.subtract(vecCoeNew, vecCoe);
		//rmse = diff.normL2();
		diff.absolute();
		rmse = diff.max();
		vecCoe = vecCoeNew;
		if(debug) 
			System.out.println("iter = " + iter + "    rmse = " + rmse);
	
		iter ++;
	  }
		return vecCoe;
	}
	
	public SimpleVector iterativeReweightedLeastSquares(int order,SimpleMatrix matX, SimpleVector vecY){
		SimpleVector vecYhat;
		SimpleVector resd, diff;
		double sigma;
		double w, rmse = 1;
		SimpleMatrix matWX = new SimpleMatrix(matX);
		SimpleVector vecWY = new SimpleVector(vecY);
		SimpleVector vecCoeNew, vecCoe;
		int maxIter = 10, iter = 0;
		
		iterativeSoftThresholding lasso = new iterativeSoftThresholding();
		double tau = 0.0005 * (1 - order/1000); 
		boolean isLasso = true;
		if(isLasso)
			 vecCoe  = lasso.runLasso(matX, vecY, tau);
		else
		     vecCoe = ordinaryLeastSquares(matX, vecY); //first estimate
		
		
	  while(iter < maxIter && rmse > 1.0e-5){
		vecYhat = SimpleOperators.multiply(matX, vecCoe);
	    resd = SimpleOperators.subtract(vecY, vecYhat);
	    sigma = madSigma(resd, matX.getCols());
	    sigma = Math.max(sigma, 1 );
		resd.divideBy(sigma);
		for(int i = 0; i < vecY.getLen(); i ++){
			w = inverseHuber(resd.getElement(i));
			//w = huber(resd.getElement(i));
			vecWY.setElementValue(i, w * vecY.getElement(i));
			for(int j = 0; j < matX.getCols(); j ++)
				matWX.setElementValue(i, j, w * matX.getElement(i, j));
		}
	
		if(isLasso)
			 vecCoeNew  = lasso.runLasso(matX, vecY, tau);
		else
		     vecCoeNew = ordinaryLeastSquares(matX, vecY);
		
		diff = SimpleOperators.subtract(vecCoeNew, vecCoe);
		
		diff.absolute();
		rmse = diff.max();
		vecCoe = vecCoeNew;

		iter ++;
	  }
		return vecCoe;
	}
	
	double getSigma(SimpleVector resd){
		double mean = 0;
		for(int i = 0; i < resd.getLen(); i++)
			mean = mean + resd.getElement(i);		
		mean = mean/resd.getLen();	
		double sigma =  1.253 * mean;
		return sigma;
	}
	
	private SimpleVector leverage(SimpleMatrix matX){
		SimpleVector h = new SimpleVector(matX.getRows());
		DecompositionQR qr = new DecompositionQR(matX);
		SimpleMatrix Q = qr.getQ();
        double temp;
		for (int i = 0; i < Q.getRows(); i ++){
			temp = Q.getRow(i).normL2();
			temp  = temp * temp;
			if (temp > 0.9999)
				h.setElementValue(i, 0.9999);
			else
				h.setElementValue(i, temp);
		}
		
		return h;
	}
	
	private SimpleVector leverageFactor(SimpleMatrix matX){
		SimpleVector h = leverage(matX);
		//System.out.println(h.toString());
		SimpleVector hfac = new SimpleVector(h.getLen());
		
		for(int i = 0; i < h.getLen(); i ++)
			hfac.setElementValue(i, 1/Math.sqrt(1 - h.getElement(i)));
		
		return hfac;
			
	}
	
	/**
	 * 
	 * @param resd: residuals
	 * @param p: the rank of matrix X, or size of X is n * p
	 * @return
	 */
	double madSigma(SimpleVector resd, int p){
		double sigma = 0;
		SimpleVector vecs = sort(resd);
		int len = resd.getLen();
		double mad; //median absolute deviation
		if((len - p)%2 == 0)
		     mad = (vecs.getElement((len + p)/2 -1) + vecs.getElement((len + p)/2 + 1))/2;
		else
			mad = vecs.getElement((len + p -1)/2);
		
		sigma = mad/0.6745;
		return sigma;
	}
	
	SimpleVector sort(SimpleVector vec){
		int len = vec.getLen();
		
		vec.absolute();
		double temp;
		for(int i = 0; i < len - 1; i ++)
			for(int j = i + 1; j < len; j ++)
				if(vec.getElement(i) > vec.getElement(j)){
					temp = vec.getElement(i);
					vec.setElementValue(i, vec.getElement(j));
					vec.setElementValue(j, temp);
				}		
		return vec;
	}
	
	/*
	 * weight functions below
	 */
	
	private double andrews(double r){
		double w = (Math.abs(r) < Math.PI ? 1 : 0) * Math.sin(r)/r;
		return w;
	}
	
	
	private double bisquare(double r){
		double w =(Math.abs(r) < 1 ? 1 : 0) * Math.pow((1 - r * r), 2);
		return w;
	}
	
	
	private double huber(double r){
		double w = 1./Math.max(1, Math.abs(r));
		return w;
	}
	
	private double fair(double r){
		double w = 1./(1 + Math.abs(r));
		return w;
	}
	
	private double cauchy(double r){
		double w = 1./(1 + r * r);
		return w;
	}
	
	private double logistic(double r){
		double w = Math.tanh(r)/r;
		return w;
	}
	
	private double inverseHuber(double r)
	{
		double w = Math.max(1, Math.abs(10*r));
		return w;
	}
	
	private double inverseFair(double r){
		double w = (1 + r * r);
		return w;
	}

}
