package VGM;

public class Spherical implements Model{
	
	double  dist;
	double sill;
	double range;
	double nug;
	
	
	public Spherical (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {
  
        double result = 0;

            double hr;
            hr = dist  / (range);
            if (dist < range) {
                result = nug + sill * hr * (1.5 - 0.5 * hr * hr);
            } else if (dist >= range) {
                result = sill+nug;
            }else if (dist==0){
            	result=0;
            }
            	
        
        return result;
	}

}
