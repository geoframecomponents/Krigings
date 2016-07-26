package VGM;

public class Periodic implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Periodic (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
		double[] result= new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                result[i] = nug + sill * (1.0 - Math.cos(2.0 * Math.PI * dist[i] / (range)));
            }
            //System.out.println(func[i]);
        }
        return result;

	}

}
