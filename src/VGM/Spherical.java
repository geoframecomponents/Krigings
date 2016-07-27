package VGM;

public class Spherical implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Spherical (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
        int length = dist.length;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            double hr;
            hr = dist[i] / (range);
            if (dist[i] < range) {
                result[i] = nug + sill * hr * (1.5 - 0.5 * hr * hr);
            } else if (dist[i] >= range) {
                result[i] = sill+nug;
            }else if (dist[i]==0){
            	result[i]=0;
            }
        }
        return result;
	}

}
