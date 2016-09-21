package theoreticalVariogram;

public class Logaritmic implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Logaritmic (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {
		double result=0;
            if (dist!= 0.0) {
                result = nug + sill * (Math.log(dist / range));
            //System.out.println(result[i]);
        }
        return result;
	}

}
