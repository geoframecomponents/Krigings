package theoreticalVariogram;

public class Periodic implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Periodic (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double  result() {
		double result= 0;

            if (dist != 0.0) {
                result  = nug + sill * (1.0 - Math.cos(2.0 * Math.PI * dist / (range)));
            }

        return result;

	}

}
