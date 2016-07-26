package VGM;

public class Power implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Power (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
        double[] result = new double[dist.length];
        for (int i = 0; i < dist.length; i++) {
            if (dist[i] != 0.0) {
                result[i] = nug + sill * (Math.pow(dist[i], range));
            }
            //System.out.println(result[i]);
        }
        return result;
	}

}
