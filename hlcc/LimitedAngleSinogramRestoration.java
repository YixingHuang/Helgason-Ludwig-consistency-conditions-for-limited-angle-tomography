package edu.stanford.rsl.tutorial.hlcc;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;
import edu.stanford.rsl.conrad.filtering.BilateralFilteringTool;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import edu.stanford.rsl.tutorial.phantoms.SheppLogan;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import java.nio.file.Paths;
/**
 * Here we perform the parallel-beam limited angle sinogram restoration based on Helgason-Ludwig consistency condition
 * @author Yixing Huang
 *
 */
public class LimitedAngleSinogramRestoration {

	public static void main(String[] args) throws Exception{
		new ImageJ();
		boolean useReconstructed = false; // if true, you can save intermediate HLCC reconstructed images for further postprocessing
		String pathReconHLCCAnalytic = "D:\\CONRAD\\CONRAD_new2\\CONRAD\\src\\edu\\stanford\\rsl\\tutorial\\hlcc\\reconHLCC.tif"; //the image of intermediate HLCC reconstruction
		boolean isRidge = false; // If ture, ridge regression is used; otherwise Lasso regression is used as default.
		boolean isPoisson = true;  //if ture, poisson noise is simulated, assuming 10e4 photons per ray.
		boolean isAnalytic = true; //compute the sinogram of shepp-logan phantom analytically instead of numerically
		boolean isShorter = false; //120 degree  if it is shorter, 120 degree is chosen instead of 160 degree limited angle tomography

		ChebyshevFourierTransformFunctions  obj = new ChebyshevFourierTransformFunctions();
		ParallelSinogramAndReconstruction  objSino = new ParallelSinogramAndReconstruction();
		
//		if(isClinical){
//		    String yourOwnDataPath = "yourOwnImag.tif";
//			ImagePlus imp0 =IJ.openImage(String yourOwnDataPath);
//			Grid2D img0 = ImageUtil.wrapImagePlus(imp0).getSubGrid(0);
//			objSino.phan = img0;
//			objSino.phan.clone().show("your own data");
//		}
		System.out.println("Many parts take time. Please wait with patience...");
		
		int sizeX = 512;
		int sizeY = 512;
		Grid2DComplex complexPhan = new Grid2DComplex(objSino.phan, false);
		complexPhan.transformForward();
		complexPhan.fftshift();
		Grid2D frequencyMag  = complexPhan.getMagnSubGrid(0, 0, sizeX, sizeY);
		frequencyMag.getGridOperator().log(frequencyMag);
		frequencyMag.show("Fourier GT log magnitude");
		
		
		objSino.maxS = objSino.imgSizeX*1.5; //
		objSino.deltaTheta = 0.5 * Math.PI / 180.0;
		objSino.maxThetaLimited = 160 * Math.PI/180;
		objSino.deltaS = 0.5;
		objSino.getCompleteSinogram();
		
		
		Grid2D sinogramComplete = objSino.sinogramComplete;
		
		if(isAnalytic){
			System.out.println("Computing analytic Shepp-Logan sinogram. It takes some time. Please wait with patience.");
			SheppLogan shep = new SheppLogan(objSino.imgSizeX, false);
			shep.setSpacing(1, 1);
			sinogramComplete = shep.analyticSinogram(objSino.maxTheta, objSino.deltaTheta, objSino.maxS, objSino.deltaS, false);
			objSino.sinogramComplete = new Grid2D(sinogramComplete);
		}
		
		Grid2D reconFull = objSino.FBPReconstructionParallel(objSino.sinogramComplete);
		//reconLimited.getGridOperator().removeNegative(reconLimited);
		reconFull.getGridOperator().removeNegative(reconFull);
		reconFull.clone().show("full sinogram reconstruction, FBP");
		
		if (isPoisson){
			System.out.print("Simulating Poisson noise with CPU. It takes some time. Please wait with patience.");
			objSino.addPoissonNoise();
			sinogramComplete = objSino.sinogramComplete;
		}
		
		sinogramComplete.clone().show("sinogram");
		
		Grid2D sinogramLimited = objSino.getLimitedAngleSinogram();
		sinogramLimited.clone().show("limited angle sinogram");
		
		Grid2D sino = new Grid2D(sinogramLimited.getSize()[0], (int)(Math.PI/objSino.deltaTheta));
		for(int i = 0; i < sino.getWidth(); i++)
			for(int j = 0; j < objSino.maxThetaLimited/objSino.deltaTheta; j++)
				sino.setAtIndex(i, j, objSino.sinogramComplete.getAtIndex(i, j));
		sino.clone().show(" limited angle sinogram");
		
		Grid2D reconLimited = objSino.FBPReconstructionParallel(sinogramLimited);
		//reconLimited.getGridOperator().removeNegative(reconLimited);
		reconLimited.clone().show("limited angle reconstruction, FBP");
	 
	
		
		
		
		Grid2D diffFBP = new Grid2D(reconFull);
			diffFBP.getGridOperator().subtractBy(diffFBP, objSino.phan);
			diffFBP.getGridOperator().abs(diffFBP);
			diffFBP.show("f FBP diff");
		/*
		Grid2D sheppMask = getSheppMask(reconLimited.getSize());
		 sheppMask.show("sheppMask");
		reconFull.getGridOperator().multiplyBy(reconFull, sheppMask);
		reconFull.clone().show("full sinogram reconstruction, FBP after mask");
		*/
		long t_start=System.currentTimeMillis();
		Grid2D reconHLCC = null;
		if(!useReconstructed){
			int maxOrder = 720;

			Grid2D momentCurveMatrixGT = new Grid2D(sinogramComplete.getSize()[1], maxOrder);
			Grid2D momentCurveMatrixRecovered = new Grid2D(sinogramComplete.getSize()[1], maxOrder);
			Grid2D sinogramRecover = new Grid2D(objSino.sinogramComplete);
			sinogramRecover.setSpacing(objSino.sinogramComplete.getSpacing());
			sinogramRecover.getGridOperator().fill(sinogramRecover, 0);

			int numThetaLimited = (int) (objSino.maxThetaLimited / objSino.deltaTheta);

			Grid1D momentCurve;
			System.out.print(
					"Start to fit the HLCC components of different orders: odd orders first and even orders afterwards!");
			System.out.print(
					"This will take some time. Please wait!");

			int n0 = 0;

			momentCurve = obj.ChebyshevFourierTransformMoment(sinogramComplete, 0);
			Grid1D momentZeroMeasured = momentCurve.getSubGrid(0, numThetaLimited);
			double mean = momentZeroMeasured.getGridOperator().sum(momentZeroMeasured) / numThetaLimited;
			System.out.println(mean);
			for (int j = 0; j < momentCurve.getNumberOfElements(); j++) {
				momentCurveMatrixGT.setAtIndex(j, n0, momentCurve.getAtIndex(j));
				momentCurveMatrixRecovered.setAtIndex(j, n0, (float) (mean));
			}
			obj.inverseChebyshevFourierTransformFromMomentCurves(momentCurveMatrixRecovered.getSubGrid(0), n0,
					sinogramRecover);

			SimpleMatrix Xodd = null, Xeven = null;
			SimpleVector momentCurveRecovered, Y;

			// odd cases
			for (int n = 1; n < maxOrder; n = n + 2) {
				System.out.print(n + ", ");
				if (n == 1)
					Xodd = obj.getMatrixXOdd(n, sinogramComplete);
				else
					Xodd = obj.getMatrixX(n, sinogramComplete, Xodd);

				momentCurve = obj.ChebyshevFourierTransformMoment(sinogramComplete, n);
				Y = new SimpleVector(momentCurve.getBuffer());

				if (isRidge) {
					float tau = 0.0005f * (1 - n / 1000);
					// momentCurveRecovered = obj.getEstimatedMomentCurveWLS(n,numThetaLimited,
					// Xodd, Y);
					momentCurveRecovered = obj.getEstimatedMomentCurveRidge(tau, numThetaLimited, Xodd, Y);
				} else // Lasso by default
					momentCurveRecovered = obj.getEstimatedMomentCurve(n, numThetaLimited, Xodd, Y);

				for (int j = 0; j < momentCurve.getNumberOfElements(); j++) {
					momentCurveMatrixGT.setAtIndex(j, n, momentCurve.getAtIndex(j));
					momentCurveMatrixRecovered.setAtIndex(j, n, (float) (momentCurveRecovered.getElement(j)));
				}
				obj.inverseChebyshevFourierTransformFromMomentCurves(momentCurveMatrixRecovered.getSubGrid(n), n,
						sinogramRecover);

			}

			System.out.println(" ");
			/**
			 * even cases
			 */
			for (int n = 2; n < maxOrder; n = n + 2) {
				System.out.print(n + ", ");
				if (n == 2)
					Xeven = obj.getMatrixXEven(n, sinogramComplete);
				else
					Xeven = obj.getMatrixX(n, sinogramComplete, Xeven);

				momentCurve = obj.ChebyshevFourierTransformMoment(sinogramComplete, n);

				Y = new SimpleVector(momentCurve.getBuffer());

				if (isRidge) {
					float tau = 0.0005f * (1 - n / 1000);
					momentCurveRecovered = obj.getEstimatedMomentCurveRidge(tau, numThetaLimited, Xeven, Y);
				} else
					momentCurveRecovered = obj.getEstimatedMomentCurve(n, numThetaLimited, Xeven, Y);

				for (int j = 0; j < momentCurve.getNumberOfElements(); j++) {
					momentCurveMatrixGT.setAtIndex(j, n, momentCurve.getAtIndex(j));
					momentCurveMatrixRecovered.setAtIndex(j, n, (float) (momentCurveRecovered.getElement(j)));
				}
				obj.inverseChebyshevFourierTransformFromMomentCurves(momentCurveMatrixRecovered.getSubGrid(n), n,
						sinogramRecover);

			}

			sinogramRecover.clone().show("sinogram Recover");

			Grid2D sinogramRecover180 = new Grid2D(sinogramRecover.getWidth(), sinogramRecover.getHeight() / 2);
			Grid2D sinoDiff = new Grid2D(sinogramRecover.getWidth(), sinogramRecover.getHeight() / 2);
			for (int i = 0; i < sinogramRecover180.getWidth(); i++) {
				for (int j = 0; j < sinogramRecover180.getHeight(); j++) {
					sinogramRecover180.setAtIndex(i, j, sinogramRecover.getAtIndex(i, j));
					sinoDiff.setAtIndex(i, j, sinogramRecover.getAtIndex(i, j) - sinogramComplete.getAtIndex(i, j));
				}
			}

			sinogramRecover180.clone().show("180-degree sinogram recovered by HLCC");
			sinoDiff.clone().show("sinogramDiff180");

			// calculate linear correlation coefficients
			Grid1D correlation = obj.getCorrCoeBetweenEstimationAndGroundTruth(momentCurveMatrixGT,
					momentCurveMatrixRecovered, objSino);
			System.out.println(" ");
			System.out.println("linear correlation coefficients: ");
			for (int i = 0; i < maxOrder; i++)
				System.out.print(correlation.getAtIndex(i) + " ");

			Grid1D sparsity = obj.getSparsityOfMomentCurves(momentCurveMatrixGT);
			System.out.println(" ");
			System.out.println("sparsity (number of coefficients smaller than epsilon): ");
			for (int i = 0; i < maxOrder; i++)
				System.out.print(sparsity.getAtIndex(i) + " ");
			System.out.println(" ");

			ImagePlus imp1 = ImageUtil.wrapGrid(momentCurveMatrixGT, null);
			IJ.saveAs(imp1, "Tiff", ("momentCurveMatrixGT.tif"));
			ImagePlus imp2 = ImageUtil.wrapGrid(momentCurveMatrixRecovered, null);
			IJ.saveAs(imp2, "Tiff", ("momentCurveMatrixRecovered.tif"));

			reconHLCC = objSino.FBPReconstructionParallel(sinogramRecover);
			reconHLCC.getGridOperator().removeNegative(reconHLCC);
			
		}
		else
		{
			ImagePlus imp = null;
			
			if(isAnalytic){	
//				if(isShorter)
//					imp = IJ.openImage(pathReconHLCCAnalytic120);
//				else if(isRidge)
//					imp = IJ.openImage(pathReconHLCCAnalyticRidge);
//				else
					imp = IJ.openImage(pathReconHLCCAnalytic);
			}
	
			reconHLCC = ImageUtil.wrapImagePlus(imp).getSubGrid(0);
		
		}
		reconHLCC.clone().show("recon with HLCC");
		
		
		System.out.println("The multiple bilateral filtering takes some time. Please wait with patience.");
		/**
		 * image fusion at the frequency domain
		 */
		boolean useCircle = true;
		imageFusionFrequencyDomain fusion = new imageFusionFrequencyDomain();
		fusion.endAngle = (int)(180 * objSino.maxThetaLimited/Math.PI);
		Grid2D imgFused = fusion.getFusedImage(reconLimited, reconHLCC, useCircle);
		imgFused.getGridOperator().removeNegative(imgFused);
		imgFused.clone().show("fused image");
		
		long t_end=System.currentTimeMillis();
		System.out.println("time is "+(t_end-t_start)/1000.0);
		
		/**
		 * image fusion, reconHLCC is filtered by a strong Bilateral filter
		 */
		
		Grid2D reconBF = new Grid2D(reconHLCC);
		//Grid2D reconBF = new Grid2D(imgFused);
		 BilateralFilteringTool bf = new  BilateralFilteringTool();
		 bf.width = 40;
		 bf.sigma_d = 30;
		 bf.sigma_r = 0.05;
		reconBF = bf.applyToolToImage(reconBF);
		reconBF = bf.applyToolToImage(reconBF);
		reconBF = bf.applyToolToImage(reconBF);
		 reconBF.show("reconstruction of HLCC after Bilateral Filter");
		
		 
		 
		 Grid2D imgFusedWithBF = fusion.getFusedImage(reconLimited, reconBF, useCircle);
		 imgFusedWithBF.getGridOperator().removeNegative(imgFusedWithBF);
		 imgFusedWithBF.clone().show("fused image with Bilateral Filter");
		 
		 /**
		  * Error measurement with respect to the ground truth phantom
		  */
		 double errLimited = obj.rmse(reconLimited, objSino.phan);
		 double errFull = obj.rmse(reconFull, objSino.phan);
		 double errHLCC = obj.rmse(reconHLCC, objSino.phan);
		 double errFused = obj.rmse(imgFused, objSino.phan);
		 double errBF = obj.rmse(reconBF, objSino.phan);
		 double errFusedBF = obj.rmse(imgFusedWithBF, objSino.phan);
		
			
		 System.out.println(" ");
		 System.out.print("errFull = " + errFull + "errLimited = " + errLimited+ " errHLCC = "+errHLCC + " errFused = " + errFused + " errBF = "+errBF+ " errFusedBF = "+errFusedBF);
	
		 
		 /**
		  * Error measurement with respect to the full sinogram reconstruction
		  */
		 
		 /**
		  * Error measurement with respect to the ground truth phantom
		  */
		 double errLimited2 = obj.rmse(reconLimited, reconFull);
		 double errFull2 = obj.rmse(reconFull, reconFull);
		 double errHLCC2 = obj.rmse(reconHLCC, reconFull);
		 double errFused2 = obj.rmse(imgFused, reconFull);
		 double errBF2 = obj.rmse(reconBF, reconFull);
		 double errFusedBF2 = obj.rmse(imgFusedWithBF, reconFull);
		
			
		 System.out.println(" ");
		 System.out.print("errFull2 = " + errFull2 + "errLimited2 = " + errLimited2+ " errHLCC2 = "+errHLCC2 + " errFused2 = " + errFused2 + " errBF2 = "+errBF2+ " errFusedBF2 = "+errFusedBF2);
	
		 
		 Grid2D diffLimited = new Grid2D(reconLimited);
			diffLimited.getGridOperator().subtractBy(diffLimited, reconFull);
			diffLimited.getGridOperator().abs(diffLimited);
			diffLimited.show("fLimitedDiff");
		 
		 Grid2D diffHLCC = new Grid2D(reconHLCC);
			diffHLCC.getGridOperator().subtractBy(diffHLCC, reconFull);
			diffHLCC.getGridOperator().abs(diffHLCC);
			diffHLCC.show("fHLCCDiff");
		 
		 Grid2D diffFused = new Grid2D(imgFused);
			diffFused.getGridOperator().subtractBy(diffFused, reconFull);
			diffFused.getGridOperator().abs(diffFused);
			diffFused.show("ffusedDiff");
		 
		 Grid2D diffBF = new Grid2D(reconBF);
		 diffBF.getGridOperator().subtractBy(diffBF, reconFull);
		 diffBF.getGridOperator().abs(diffBF);
		 diffBF.show("fBFDiff");
		 
		 Grid2D diffFused2 = new Grid2D(imgFusedWithBF);
			diffFused2.getGridOperator().subtractBy(diffFused2, reconFull);
			diffFused2.getGridOperator().abs(diffFused2);
			diffFused2.show("ffused2Diff");
		 
	}
	
}
