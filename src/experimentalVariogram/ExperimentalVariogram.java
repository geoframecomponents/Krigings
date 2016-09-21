/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
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
package experimentalVariogram;

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
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

@Description("Experimental semivariogram algorithm.")
@Documentation("Experimental semivariogram")
@Author(name = "Giuseppe Formetta, Francesco Adami, Silvia Franceschi & Marialaura Bancheri")
@Keywords("Experimental semivariogram, Kriging, Hydrology")
@Label(JGTConstants.STATISTICS)
@Name("variogram")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class ExperimentalVariogram extends JGTModel {


	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;


	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;


	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;


	@Description("The file with the measured data, to be interpolated.")
	@In
	public HashMap<Integer, double[]> inData = null;

	
	@Description("Spatial separation distance up to which point pairs are included in semivariance estimates; "
			+ "as a default, the length of the diagonal of the box spanning the data is divided by three.")
	@In
	public double pCutoff;

	
	@Description("The Experimental Distances.")
	@Out
	public HashMap<Integer, double[]>  outDistances;


	@Description("The Experimental Variogram.")
	@Out
	public HashMap<Integer, double[]>  outExperimentalVariogram;
	

	boolean areAllEquals;
	

	int differents;


	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	private HortonMessageHandler msg = HortonMessageHandler.getInstance();


	/**
	 * Process.
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {

		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();

		/*
		 * counter for the number of station with measured value !=0.
		 */
		int n1 = 0;

		/*
		 * Store the station coordinates and measured data in the array.
		 */
		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid))
						.intValue();
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
				Coordinate coordinate = ((Geometry) feature
						.getDefaultGeometry()).getCentroid().getCoordinate();


				// h is the vector with measured data
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {

					/*
					 * skip data for non existing stations.
					 * Also skip novalues.
					 */
					continue;
				}

				if (Math.abs(h[0]) >= 0.0) { // TOLL
					xStationList.add(coordinate.x);
					yStationList.add(coordinate.y);
					zStationList.add(z);
					hStationList.add(h[0]);
					n1 = n1 + 1;
				}

			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * The coordinates of the station points plus in the last position a place
		 * for the coordinate of the point to interpolate.
		 */
		double[] xStation = new double[nStaz];
		double[] yStation = new double[nStaz];
		double[] zStation = new double[nStaz];
		double[] hStation = new double[nStaz];

		areAllEquals = true;

		differents = 0;

		if (nStaz != 0) {
			xStation[0] = xStationList.get(0);
			yStation[0] = yStationList.get(0);
			zStation[0] = zStationList.get(0);
			hStation[0] = hStationList.get(0);

			double previousValue = hStation[0];

			for (int i = 1; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation(
						xStation, yStation, zStation, hStation, xTmp, yTmp,
						zTmp, hTmp, i, false, pm);

				if (!doubleStation) {
					xStation[i] = xTmp;
					yStation[i] = yTmp;
					zStation[i] = zTmp;
					hStation[i] = hTmp;
					if (areAllEquals && hStation[i] != previousValue) {

						areAllEquals = false;
					}
					if (hStation[i] != previousValue) {
						differents += 1;
					}
					previousValue = hStation[i];
				}
			}
		}


		// number of different stations
		if (differents > 2) {
			

			double[][] outResult = processAlgorithm(xStation, yStation, hStation, pCutoff);
			storeResult(outResult);

		}
		else{
			System.out.println("Only 1 data >0 or All the data are equal. Variogram is not running");
		}

	}



	/**
	 * Process algorithm.
	 *
	 * @param xStation the vector containing the x value of the station
	 * @param yStation the vector containing the y value of the station
	 * @param hStation the vector containing the variable value of the station
	 * @param Cutoffinput the cutoff input
	 * @return the double[][] matrix of the results of the processing
	 */
	public double[][] processAlgorithm(double[] xStation, double yStation[],
			double[] hStation, double Cutoffinput) {

		double x1, x2, y1, y2;
		double dDifX, dDifY;
		double value;
		double mean = 0;
		double maxDistance = 0;


		double Cutoff;
		int iCount = xStation.length;
		double distanceMatrix[][] = new double[iCount][iCount];
		double x_max = xStation[0], y_max = yStation[0], diagonal;
		double x_min = xStation[0], y_min = yStation[0];
		for (int i = 1; i < iCount; i++) {

			x_min = Math.min(x_min, xStation[i]);
			y_min = Math.min(y_min, yStation[i]);
			x_max = Math.max(x_max, xStation[i]);
			y_max = Math.max(y_max, yStation[i]);

		}

		diagonal = Math.sqrt((x_max - x_min) * (x_max - x_min)
				+ (y_max - y_min) * (y_max - y_min));

		if (Cutoffinput == 0) {
			Cutoff = diagonal / 3;
		} else
			Cutoff = Cutoffinput;


		// Compute the distance matrix
		for (int i = 0; i < iCount; i++) {
			x1 = xStation[i];
			y1 = yStation[i];
			value = hStation[i];

			mean += value;

			for (int j = 0; j < iCount; j++) {

				x2 = xStation[j];
				y2 = yStation[j];

				dDifX = x2 - x1;
				dDifY = y2 - y1;

				// Pitagora theorem
				distanceMatrix[i][j] = Math.sqrt(dDifX * dDifX + dDifY * dDifY); 

				maxDistance = Math.max(maxDistance, distanceMatrix[i][j]);

			}
		}

		// compute the mean of the input hStation
		mean /= (double) iCount; 

		double[][] result = calculate(Cutoff, distanceMatrix,hStation,  mean, maxDistance);
		
		return result;

	} 

	/**
	 * Calculate the variances and the distances
	 *
	 * @param cutoff the cutoff
	 * @param distanceMatrix the distance matrix
	 * @param the vector containing the variable value of the station
	 * @param mean the mean value of the input data
	 * @param maxDistance the max distance value
	 * @return the double[][] matrix with the results (the variances and the distances)
	 */
	public double[][] calculate( double cutoff, double[][] distanceMatrix, double[] hStation, double mean,
			double maxDistance) {


		int Cutoff_divide = 15;
		double binAmplitude = cutoff / Cutoff_divide;


		// number of distance for each bin
		int iClasses = (int) (maxDistance / binAmplitude + 2); 

		// definition of the vectors containing the variance, covariance, semivariance, 
		// 	number of the points in the specified bin..


		double[] m_dSemivar = new double[iClasses]; 

		int[] iPointsInClass = new int[iClasses]; 


		double[] m_ddist = new double[iClasses];
		int NonZero=0;

		for (int i = 0; i < distanceMatrix.length; i++) {
			
			// first cycle input hStation 
			double value1 = hStation[i]; 

			for (int j = i + 1; j < distanceMatrix.length; j++) {
				if (distanceMatrix[i][j] > 0 && distanceMatrix[i][j] < cutoff) {

					// return the class of considered distance 
					int iClass = (int) Math.floor((distanceMatrix[i][j]) / binAmplitude);
					
					// counts the number of distances of each class
					iPointsInClass[iClass]++; 

					// second cycle input hStation
					double value2 = hStation[j]; 

					// compute the numerator of the semivariance 
					double dSemivar = Math.pow((value1 - value2), 2.); 

					// sum all the semivariances for the considered class
					
					NonZero=(m_dSemivar[iClass]!=0)?NonZero:NonZero+1;
					m_dSemivar[iClass] += dSemivar; 

					m_ddist[iClass] += distanceMatrix[i][j]; 
					

				}
			}


		}

		double[][] result = new double[NonZero][2];

		for (int i = 0; i < NonZero; i++) {
				// Compute the semivariance
				m_dSemivar[i] /= (2. * iPointsInClass[i]); 

				// Compute the mean distance for each class
				m_ddist[i] /= iPointsInClass[i];

				result[i][0] = m_ddist[i];
				result[i][1] = m_dSemivar[i];	
		}
	
		return result;

	}
	
	/**
	 * Store result.
	 *
	 * @param result are the resulting variances and the distances
	 * @throws SchemaException the schema exception
	 */
	private void storeResult(double [][] result) 
			throws SchemaException {
		outDistances = new HashMap<Integer, double[]>();
		outExperimentalVariogram = new HashMap<Integer, double[]>();
		
		for (int i=0;i<result.length;i++){
		outDistances.put(i, new double[]{result[i][0]});
		outExperimentalVariogram.put(i, new double[]{result[i][1]});
		}
	}

}