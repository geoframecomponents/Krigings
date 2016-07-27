package krigings;

public class Spline implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Spline (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {
        double  result = 0;

            if (dist < range) {
                result = nug + sill * (dist * dist * Math.log(dist));
            }else if (dist>= range) {
                result = sill+nug;
            }

        return result;
	}

}
