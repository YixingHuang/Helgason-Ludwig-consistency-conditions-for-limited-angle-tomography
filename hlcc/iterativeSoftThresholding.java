package edu.stanford.rsl.tutorial.hlcc;
import edu.stanford.rsl.conrad.numerics.DecompositionSVD;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix.MatrixNormType;

/**
 * Iterative thresholding algorithm, 
 * I. Daubechies, M. Defrise, and C. De Mol. 
 * An iterative thresholding algorithm for linear inverse problems with a sparsity constraint.
 * Communications on pure and applied mathematics, 57(11):1413–1457, 2004.
 * From Xiaolin Huang's Matlab code
 * @author Yixing Huang
 *
 */
public class iterativeSoftThresholding {
	public double eps = 1.0e-4;
	public iterativeSoftThresholding(){
		
	}
	
	public SimpleVector runLasso(SimpleMatrix Phi,SimpleVector Y,  double tau){
        double s = 1.5;
		//double Xnorm = X.norm(MatrixNormType.MAT_NORM_L2);
		DecompositionSVD decompositionSVD = new DecompositionSVD(Phi);
		double Xnorm = decompositionSVD.getS().getElement(0, 0);
		Phi.divideBy(Xnorm*s);

		Y.divideBy(Xnorm*s);
		
		SimpleMatrix PhiT =Phi.transposed();
		SimpleVector x_k_1 = SimpleOperators.multiply(PhiT, Y);
		SimpleVector x_k = new SimpleVector(x_k_1.getLen());
		SimpleVector residual;
		double diff = SimpleOperators.subtract(x_k_1, x_k).normL2()/x_k.normL2();
		while(diff>eps){
			x_k = (SimpleVector) x_k_1.clone();
			residual = SimpleOperators.subtract(Y, SimpleOperators.multiply(Phi, x_k));
			x_k_1 = SimpleOperators.add(x_k, SimpleOperators.multiply(PhiT, residual));
			x_k_1 = soft_threshold(x_k_1,tau);
			diff = SimpleOperators.subtract(x_k_1, x_k).normL2()/x_k.normL2();
		}
	
		return x_k_1;
	}
	
	public SimpleVector runWithInitialization(SimpleMatrix Phi, SimpleVector Y, double tau, SimpleVector xini){
		double s = 1.5;
		DecompositionSVD decompositionSVD = new DecompositionSVD(Phi);
		double Xnorm = decompositionSVD.getS().getElement(0, 0);
		Phi.divideBy(Xnorm * s);
		Y.divideBy(Xnorm *s);
		
		SimpleMatrix PhiT =Phi.transposed();
		SimpleVector x_k_1 = SimpleOperators.multiply(PhiT, Y);
		SimpleVector x_k = new SimpleVector(x_k_1.getLen());
		SimpleVector residual;
		x_k_1 = (SimpleVector) xini.clone();
		double diff = SimpleOperators.subtract(x_k_1, x_k).normL2()/x_k.normL2();
		while(diff>eps){
			x_k = (SimpleVector) x_k_1.clone();
			residual = SimpleOperators.subtract(Y, SimpleOperators.multiply(Phi, x_k));
			x_k_1 = SimpleOperators.add(x_k, SimpleOperators.multiply(PhiT, residual));
			x_k_1 = soft_threshold(x_k_1,tau);
			diff = SimpleOperators.subtract(x_k_1, x_k).normL2()/x_k.normL2();
		}
	
		return x_k_1;
	}
	private SimpleVector soft_threshold(SimpleVector x, double tau){
		double val;
		SimpleVector y = new SimpleVector(x.getLen());
		for(int i = 0; i < x.getLen(); i++){
			val = x.getElement(i);
			if (Math.abs(val)<tau)
				y.setElementValue(i, 0);
			else if(val<-tau)
				y.setElementValue(i, val + tau);
			else
				y.setElementValue(i, val - tau);
		
		}
		return y;
	}

}
