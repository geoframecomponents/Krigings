/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package trendAnalysis;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;
import org.opengis.feature.simple.SimpleFeature;


import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import flanagan.analysis.Regression;

@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Marialaura Bancheri, Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello")
@Keywords("Kriging, Hydrology")
@Label("")
@Name("kriging")
@Status()
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class TrendAnalysis extends JGTModel {



	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;


	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;


	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;


	@Description("The file with the measured data to be interpolated.")
	@In
	public HashMap<Integer, double[]> inData = null;


	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();


	@Description("Include zeros in computations (default is true).")
	@In
	public boolean doIncludezero = true;


	@Description("The threshold on correlation coefficient for the trend in detrendend mode.")
	@In
	public double thresholdCorrelation;

	@Description("Degree of polynomial regression, default is 1")
	@In
	public int regressionOrder=1;


	@Description("Switch for detrended mode.")
	@In
	@Out
	public boolean doDetrended = true;



	@Description("Intercept of the linear regression")
	@In
	@Out
	public double trend_intercept=0;

	@Description("Coefficient of the linear regression")
	@In
	@Out
	public double trend_coefficient=0;


	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outResiduals= null;


	private HortonMessageHandler msg = HortonMessageHandler.getInstance();



	int id;



	@Execute
	public void execute() throws Exception {

		verifyInput();

		// create the arraylist containing the station with the measurements 
		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();
		List<Integer> idStationList = new ArrayList<Integer>();

		/*
		 * counter for the number of station with measured value !=0.
		 */
		int n1 = 0;


		/*
		 * Store the station coordinates and measured data in the array.
		 * Skip data for non existing stations and also skip novalues.
		 */

		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid)).intValue();

				double z = 0;
				if (fStationsZ != null) {
					try {
						z = ((Number) feature.getAttribute(fStationsZ))
								.doubleValue();
					} catch (NullPointerException e) {
						pm.errorMessage(msg.message("kriging.noStationZ"));
						throw new Exception(msg.message("kriging.noStationZ"));

					}
				}


				Coordinate coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {

					continue;
				}

				if (doIncludezero) {
					if (Math.abs(h[0]) >= 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);
						n1 = n1 + 1;
					}

				} else {
					if (Math.abs(h[0]) > 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);

						n1 = n1 + 1;
					}
				}
			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * Check if the coordinates or the values are the same for all the measurements stations.
		 * xStationInitialSet has the dimensions of the coordinates of the measurements points 
		 * 
		 */


		double[] xStationInitialSet = new double[nStaz];
		double[] yStationInitialSet = new double[nStaz];
		double[] zStationInitialSet = new double[nStaz];
		double[] hStationInitialSet = new double[nStaz];
		int[] idStationInitialSet = new int[nStaz];

		boolean areAllEquals = true;
		if (nStaz != 0) {
			xStationInitialSet[0] = xStationList.get(0);
			yStationInitialSet[0] = yStationList.get(0);
			zStationInitialSet[0] = zStationList.get(0);
			hStationInitialSet[0] = hStationList.get(0);
			idStationInitialSet[0] = idStationList.get(0);
			double previousValue = hStationInitialSet[0];

			/* for each station added to the vector, it checks the coordinate/values and if they are
			 * not null or equal, it adds to the list of the availble stations. If the coordinates
			 * or the values are all different, the flag areAllEquals  becomes false.   
			 */

			for (int i = 0; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);
				int idTmp = idStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation( xStationInitialSet, yStationInitialSet, zStationInitialSet, 
						hStationInitialSet,xTmp,yTmp, zTmp, hTmp, i, false, pm);
				if (!doubleStation) {
					xStationInitialSet[i] = xTmp;
					yStationInitialSet[i] = yTmp;
					zStationInitialSet[i] = zTmp;
					hStationInitialSet[i] = hTmp;
					idStationInitialSet[i] = idTmp;
					if (areAllEquals && hStationInitialSet[i] != previousValue) {
						areAllEquals = false;
					}
					previousValue = hStationInitialSet[i];
				}
			}
		}


		double[] hresiduals=null;


		if (n1 != 0) {
			if (!areAllEquals && n1 > 1) {

				Regression r = new Regression();

				r = new Regression(zStationInitialSet, hStationInitialSet);
				r.polynomial(regressionOrder);


				/*If there is a trend for meteorological
				 * variables and elevation and it is statistically significant 
				 * then the residuals from this linear trend
				 * are computed for each meteorological stations.
				 */
				if (Math.abs(r.getXYcorrCoeff()) > thresholdCorrelation) {

					trend_intercept=r.getBestEstimates()[0];
					trend_coefficient=r.getBestEstimates()[1];
					hresiduals = r.getResiduals();

				} else {
					System.out.println("The trend is not significant");
					doDetrended=false;
					hresiduals=hStationInitialSet;

				}
			}


		}else{
			pm.errorMessage("No value for this time step");

		}

		storeResult(hresiduals,idStationInitialSet);

	}


	/**
	 * Verify the input of the model.
	 */
	private void verifyInput() {
		if (inData == null || inStations == null) {
			throw new NullPointerException( msg.message("kriging.stationProblem"));
		}

	}


	/**
	 * Round.
	 *
	 * @param value is the value of the variable considered 
	 * @param places the places to consider after the comma
	 * @return the double value of the variable rounded
	 */
	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		long factor = (long) Math.pow(10, places);
		value = value * factor;
		long tmp = Math.round(value);
		return (double) tmp / factor;
	}



	/**
	 * Store the result in a HashMap
	 *
	 * @param hresiduals the vector with the residuals
	 * @param hStation the vector with the values of the stations
	 * @param id the associated id of the calculating points.
	 * @throws SchemaException the schema exception
	 */
	private void storeResult(double [] hresiduals, int [] id ) throws SchemaException {
		outResiduals = new HashMap<Integer, double[]>();
		for (int i = 0; i < hresiduals.length; i++) {			
			hresiduals[i]=(hresiduals[i]<0)?0:hresiduals[i];
			outResiduals.put(id[i], new double[] { hresiduals[i] });
		}
	}




}
