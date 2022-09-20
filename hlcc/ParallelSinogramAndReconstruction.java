package edu.stanford.rsl.tutorial.hlcc;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.filtering.PoissonNoiseFilteringTool;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import edu.stanford.rsl.tutorial.parallel.ParallelBackprojector2D;
import edu.stanford.rsl.tutorial.parallel.ParallelProjector2D;
import edu.stanford.rsl.tutorial.phantoms.Phantom;
import edu.stanford.rsl.tutorial.phantoms.SheppLogan;

public class ParallelSinogramAndReconstruction {
	public int imgSizeX = 512, imgSizeY = imgSizeX;
	public ParallelProjector2D projector;
	public ParallelBackprojector2D backproj;
	
	public double deltaTheta = 0.5 * Math.PI / 180.0;
	public double maxTheta = 360 * Math.PI / 180.0;
	public double maxS = imgSizeX * 1f;
	public double deltaS = 0.5f;
	public double maxThetaLimited = 150 * Math.PI/180.0;
	public Grid2D phan;
	public Grid2D sinogramComplete;
	public double spacingX = 1, spacingY = 1;
	
	public ParallelSinogramAndReconstruction(){
		this.phan = new SheppLogan(imgSizeX, false);
		this.phan.setOrigin(imgSizeX/2, imgSizeY/2);
		this.phan.setSpacing(spacingX, spacingY);
		//this.phan.show("phantom");	
	}
	
	public void getCompleteSinogram(){
		 projector = new ParallelProjector2D(maxTheta, deltaTheta, maxS, deltaS);
		 //sinogramComplete = projector.projectRayDrivenCL(this.phan);
		 sinogramComplete = projector.projectRayDriven(this.phan);
		
	}
	
	public void addPoissonNoise() throws Exception{
		float amp = 30;
		sinogramComplete.getGridOperator().divideBy(sinogramComplete, amp);
		sinogramComplete.clone().show("sinogram");
		
		Grid2D sinogramIdeal0 = new Grid2D(sinogramComplete);
		
		float I0 = 1.e4f;
		double val;
		Grid2D I = new Grid2D(sinogramComplete.getWidth(), sinogramComplete.getHeight());
		for (int i = 0; i < sinogramComplete.getWidth(); i ++)
			for(int j = 0; j < sinogramComplete.getHeight(); j++)
			{
				val = I0 * Math.pow(Math.E, -sinogramComplete.getAtIndex(i, j));
				I.setAtIndex(i, j, (float)(val));
			}
		
		I.clone().show("I before Poisson");
		PoissonNoiseFilteringTool poisson = new PoissonNoiseFilteringTool();
		poisson.applyToolToImage(I);
		I.clone().show(" I after possion noise");
		
		for (int i = 0; i < sinogramComplete.getWidth(); i ++)
			for(int j = 0; j < sinogramComplete.getHeight(); j++)
			{
				val = - Math.log(I.getAtIndex(i, j)/I0);
				sinogramComplete.setAtIndex(i, j, (float)(val));
			
			}
	
		sinogramIdeal0.getGridOperator().subtractBy(sinogramIdeal0, sinogramComplete);
		sinogramIdeal0.clone().show("Poisson noise");
		sinogramComplete.getGridOperator().multiplyBy(sinogramComplete, amp);
	}
	
	public Grid2D getLimitedAngleSinogram(){
		int angleIndex = (int) (this.maxThetaLimited/this.deltaTheta);
		Grid2D sinogram = new Grid2D(sinogramComplete);
		sinogram.setSpacing(sinogramComplete.getSpacing());
		for(int i = 0; i < sinogram.getSize()[0]; i ++)
			for(int j = angleIndex; j < Math.PI/this.deltaTheta; j ++){
				sinogram.setAtIndex(i, j, 0);
				sinogram.setAtIndex(i, (int)(j + Math.PI/this.deltaTheta), 0);
			}
		for(int i = 0; i < sinogram.getSize()[0]; i ++)
			for(int j = 0; j < (Math.PI - this.maxThetaLimited)/this.deltaTheta; j ++){
				sinogram.setAtIndex(i, j, 0);
				sinogram.setAtIndex(i, (int)(j + Math.PI/this.deltaTheta), 0);
			}
		return sinogram;
	}
	
	public Grid2D FBPReconstructionParallel(Grid2D sinogram){
		Grid2D filteredSinogram = new Grid2D(sinogram);
		RamLakKernel ramLak = new RamLakKernel(sinogram.getWidth(), 1);
		for (int i = 0; i < sinogram.getHeight(); i ++)
			ramLak.applyToGrid(filteredSinogram.getSubGrid(i));
		filteredSinogram.setSpacing(sinogram.getSpacing()[0], sinogram.getSpacing()[1]);
		backproj = new ParallelBackprojector2D(imgSizeX, imgSizeY, (float) this.spacingX, (float) this.spacingY);
		Grid2D recon = backproj.backprojectPixelDriven(filteredSinogram);
		recon.getGridOperator().multiplyBy(recon, 2.f);
		return recon;
	}
}
