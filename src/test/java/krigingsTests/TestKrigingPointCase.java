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

package krigingsTests;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geotools.data.simple.SimpleFeatureCollection;

import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import org.junit.Test;
import static org.junit.Assert.*;

import krigingsPointCase.Krigings;


public class TestKrigingPointCase{


	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which all the station have the same value.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */
	
	@Test
	public void testKriging2() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/PointCase/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/krigings/PointCase/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="resources/Input/krigings/PointCase/rain_test2A_allNoValue.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();

		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";
        kriging.maxdist=40368.0;

        kriging.range = 123537.0;
        kriging.nugget = 0.0;
        kriging.sill= 1.678383;
        kriging.pSemivariogramType="linear";
        

		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/krigings/PointCase/kriging_interpolated_NoValue.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		//
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
	
	
			HashMap<Integer, double[]> result = kriging.outData;
			
			/*
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iterator = pointsToInterpolateResult.iterator();
			while( iterator.hasNext() ) {
				int id = iterator.next();
				double[] actual = result.get(id);
				assertEquals(1.0, actual[0], 0);
			}*/
			
			writer.inData = result;
			writer.writeNextLine();
		}
			
		//
		reader.close();
		writer.close();
	}
	// /////////////////////////////////////////////////////////////////////////////////////////
	// ///////////////////////////////FINE TEST 2
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	//
	// /////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////// TEST 3
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	// /**
	// * Run the kriging models.
	// *
	// * <p>
	// * This is the case that defaultMode=0.
	// * </p>
	// * @throws Exception
	// * @throws Exception
	// */
	public void testKriging4() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/PointCase/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/krigings/PointCase/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/PointCase/rain_test.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();
		//kriging.pm = pm;
		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";

		
		kriging.pSemivariogramType="linear";

        kriging.range = 123537.0;
        kriging.nugget = 0.0;
        kriging.sill= 1.678383;
        kriging.maxdist=1000;



		//
		kriging.doIncludezero = false;
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file =  "resources/Output/krigings/PointCase/kriging_interpolated_2.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		//
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = kriging.outData;
			

			//
			writer.inData = result;
			writer.writeNextLine();
		}
		//
		reader.close();
		writer.close();
	}



	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which there is only one station.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */
	public void testKriging5() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/PointCase/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/krigings/PointCase/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/PointCase/rain_test3A.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();
		//kriging.pm = pm;
		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";

		
		 // Set up the model in order to use the variogram with an explicit integral scale and
        //variance.
		 
		kriging.pSemivariogramType="linear";
        kriging.range = 123537.0;
        kriging.nugget = 0.0;
        kriging.sill= 1.678383;
        //kriging.maxdist=1000;

		 
		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/krigings/PointCase/kriging_interpolated_3.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		int j = 0;
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			
			 // Extract the result.
			 
			HashMap<Integer, double[]> result = kriging.outData;
			
			
			
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
			double expected;
			if (j == 0) {
				expected = 10.0;
			} else if (j == 1) {
				expected = 15;
			} else if (j == 2) {
				expected = 1;
			} else if (j == 3) {
				expected = 2;
			} else if (j == 4) {
				expected = 2;
			} else if (j == 5) {
				expected = 0;
			} else if (j == 6) {
				expected = 0;
			} else if (j == 7) {
				expected = 23;
			} else if (j == 8) {
				expected = 50;
			} else if (j == 9) {
				expected = 70;
			} else if (j == 10) {
				expected = 30;
			} else if (j == 11) {
				expected = 10;
			} else if (j == 12) {
				expected = 2;
			} else {
				expected = 1.0;
			}
			//
			while( iteratorTest.hasNext() ) {
				int id = iteratorTest.next();
				double[] actual = result.get(id);
				//
				assertEquals(expected, actual[0], 0);
			}
			
			
			writer.inData = result;
			writer.writeNextLine();
			j++;
		}
		//
		reader.close();
		writer.close();
	}

}
