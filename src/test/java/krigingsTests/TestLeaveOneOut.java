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

import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import org.junit.Test;

import leaveOneOut.Krigings;

import static org.junit.Assert.*;


public class TestLeaveOneOut{


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
		stationsReader.file = "resources/Input/krigings/PointCase/leave.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//

		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/PointCase/rain_test1.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		 reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();

		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		kriging.fStationsZ="QUOTA";		
		//

		 kriging.range = 123537.0;
		 kriging.nugget = 0.0;
		 kriging.sill= 1.678383;
		 kriging.pSemivariogramType="linear";
        
		//kriging.doDetrended=true;
		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/krigings/PointCase/kriging_interpolated_leave.csv";
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

			writer.inData = result;
			writer.writeNextLine();
		}
		//
		reader.close();
		writer.close();
	}
}