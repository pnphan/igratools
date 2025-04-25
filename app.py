from flask import Flask, render_template, jsonify, send_file, request
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc
import requests
import zipfile
from urllib.parse import urljoin
import shutil
import io

app = Flask(__name__)

def download_locations():
    """
    Download and parse the IGRA station list to get station locations.
    
    Returns:
    --------
    dict
        A dictionary with station IDs as keys and (latitude, longitude) tuples as values.
    """
    # URL for the IGRA station list
    url = 'https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/doc/igra2-station-list.txt'
    
    # Get the station list without downloading the file
    response = requests.get(url)
    response.raise_for_status()  # Raise an exception for HTTP errors
    
    # Parse the station list
    stations = []
    for line in response.text.splitlines():
        if line.strip():  # Skip empty lines
            try:
                # Station ID is in columns 1-11
                station_id = line[:11].strip()
                
                # Latitude is in columns 13-20
                latitude = float(line[12:20].strip())
                
                # Longitude is in columns 22-30
                longitude = float(line[21:30].strip())
                
                # Station name is in columns 42-71
                station_name = line[41:71].strip()
                
                # Skip stations with invalid coordinates (missing values in the file)
                if latitude != -99.9999 and longitude != -999.9999:
                    stations.append({
                        'id': station_id,
                        'name': station_name,
                        'lat': latitude,
                        'lon': longitude
                    })
            except (ValueError, IndexError):
                # Skip lines that can't be parsed
                continue
    
    return stations

def get_availability(station_id, download=True, download_dir=None):
    """
    Get the availability of IGRA data for a given station.
    """
    availability = []
    try:
        if download:
            base_url = "https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/"

            # Construct the URL for the station's zip file
            zip_filename = f"{station_id}-data.txt.zip"
            zip_url = urljoin(base_url, zip_filename)
            
            try:
                # Download the zip file directly to memory without saving it
                print(f"Downloading data for station {station_id}...")
                response = requests.get(zip_url, stream=True)
                response.raise_for_status()  # Raise an exception for HTTP errors
                
                # Create a BytesIO object to hold the zip file in memory
                zip_buffer = io.BytesIO(response.content)
                
                # Extract and process the text file directly from memory
                print(f"Processing data for station {station_id}...")
                with zipfile.ZipFile(zip_buffer) as zip_ref:
                    # Get the name of the text file in the zip
                    txt_file_name = None
                    for file_info in zip_ref.infolist():
                        if file_info.filename.endswith('.txt'):
                            txt_file_name = file_info.filename
                            break
                    
                    if txt_file_name:
                        # Read the text file directly from the zip archive
                        with zip_ref.open(txt_file_name) as file:
                            # Decode bytes to string and process line by line
                            for line in io.TextIOWrapper(file, encoding='utf-8'):
                                line = line.strip()
                                if line.startswith('#'):
                                    year = int(line[13:17])
                                    month = int(line[18:20])
                                    day = int(line[21:23])
                                    hour = int(line[24:26])
                                    availability.append([year, month, day, hour])
                    else:
                        raise Exception("No text file found in the zip archive")
                
                print(f"Successfully processed data for {station_id}")
            
            except requests.exceptions.RequestException as e:
                print(f"Error downloading data for station {station_id}: {e}")
                return None
            except zipfile.BadZipFile:
                print(f"Error: The downloaded file for {station_id} is not a valid zip file.")
                return None
            except Exception as e:
                print(f"Unexpected error processing data for {station_id}: {e}")
                return None
        else:
            # If download is False, we still need the download_dir to exist and contain the data file
            if download_dir is None:
                download_dir = os.path.join(os.getcwd(), str(datetime.now().strftime("%Y-%m-%d")))

            data_file = os.path.join(download_dir, f"{station_id}-data.txt")

            if not os.path.exists(data_file):
                print(f"Error: Data file for station {station_id} not found at {data_file}")
                return None
            
            with open(data_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('#'):
                        year = int(line[13:17])
                        month = int(line[18:20])
                        day = int(line[21:23])
                        hour = int(line[24:26])
                        availability.append([year, month, day, hour])
                    
        return availability
                    
    except Exception as e:
        print(f"Error fetching availability for station {station_id}: {e}")
        return None


def netcdf_app(station_name, download=True, download_dir=None, main=False):
    """
    Create a NetCDF file from IGRA data for a given station.
    
    Parameters:
    -----------
    station_list : list of str
        A list of IGRA station identifiers (e.g., ['USM00072401', 'USM00072421'])
    download : bool, optional
        If True, the IGRA data files will be downloaded. If False, the IGRA data files will not be downloaded.
    download_dir : str, optional
        Path to the directory to save the IGRA data files. If None, a default name will be used.
        
    Returns:
    --------
    str
        Path to the created NetCDF file, or None if an error occurred
    """
    try:
                    
        if download_dir is None:
            download_dir = os.path.join(os.getcwd(), 'igra-pytools', str(datetime.now().strftime("%Y-%m-%d")))

        # Create the directory if it doesn't exist
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

        # Create the NetCDF file
        if main:
            output_file = os.path.join(download_dir, f"{station_name}-main.nc")
        else:
            output_file = os.path.join(download_dir, f"{station_name}-full.nc")


        if download:
            base_url = "https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/"

            # Construct the URL for the station's zip file
            zip_filename = f"{station_name}-data.txt.zip"
            zip_url = urljoin(base_url, zip_filename)
            
            try:
                # Download the zip file directly to memory without saving it
                print(f"Downloading data for station {station_name}...")
                response = requests.get(zip_url, stream=True)
                response.raise_for_status()  # Raise an exception for HTTP errors
                
                # Create a BytesIO object to hold the zip file in memory
                zip_buffer = io.BytesIO(response.content)
                
                # Extract only the text file directly from memory to disk
                print(f"Extracting data for station {station_name}...")
                with zipfile.ZipFile(zip_buffer) as zip_ref:
                    # Get the name of the text file in the zip
                    txt_file_name = None
                    for file_info in zip_ref.infolist():
                        if file_info.filename.endswith('.txt'):
                            txt_file_name = file_info.filename
                            break
                    
                    if txt_file_name:
                        # Extract only the text file
                        with zip_ref.open(txt_file_name) as source, open(os.path.join(download_dir, f"{station_name}-data.txt"), 'wb') as target:
                            shutil.copyfileobj(source, target)
                    else:
                        raise Exception("No text file found in the zip archive")
                
                print(f"Successfully downloaded and extracted data for {station_name} to {download_dir}")
            
            except requests.exceptions.RequestException as e:
                print(f"Error downloading data for station {station_name}: {e}")
                return None
            except zipfile.BadZipFile:
                print(f"Error: The downloaded file for {station_name} is not a valid zip file.")
                return None
            except Exception as e:
                print(f"Unexpected error processing data for {station_name}: {e}")
                return None

        data_file = os.path.join(download_dir, f"{station_name}-data.txt")


        if not os.path.exists(data_file):
            print(f"Error: Data file for station {station_name} not found at {data_file}")
            return None


        # First pass: count soundings and find max levels
        sounding_count = 0
        max_levels = 0
        current_levels = 0
        
        with open(data_file, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#'):
                    sounding_count += 1
                    if current_levels > max_levels:
                        max_levels = current_levels
                    current_levels = 0
                else:
                    current_levels += 1
            
            # Check the last sounding
            if current_levels > max_levels:
                max_levels = current_levels
        
        print(f"Found {sounding_count} soundings with maximum {max_levels} levels")
        
        # Create arrays to store data
        dates = np.zeros(sounding_count, dtype=np.int32)
        times = np.zeros(sounding_count, dtype=np.int32)
        reltimes = []
        numlevs = np.zeros(sounding_count, dtype=np.int32)
        p_srcs = []
        np_srcs = []

        pressure = np.full((sounding_count, max_levels), np.nan)
        gph = np.full((sounding_count, max_levels), np.nan)
        temp = np.full((sounding_count, max_levels), np.nan)
        rh = np.full((sounding_count, max_levels), np.nan)
        wspd = np.full((sounding_count, max_levels), np.nan)
        wdir = np.full((sounding_count, max_levels), np.nan)
        lvltyp1 = np.full((sounding_count, max_levels), np.nan)
        lvltyp2 = np.full((sounding_count, max_levels), np.nan)
        etime = np.full((sounding_count, max_levels), np.nan)
        pflag = np.full((sounding_count, max_levels), np.nan)
        zflag = np.full((sounding_count, max_levels), np.nan)
        tflag = np.full((sounding_count, max_levels), np.nan)
        dpdp = np.full((sounding_count, max_levels), np.nan)
        
        
        # Second pass: read data
        sounding_idx = -1
        level_idx = 0
        
        with open(data_file, 'r') as file:
            for line in file:
                line = line.strip()
                
                if line.startswith('#'):
                    sounding_idx += 1
                    level_idx = 0
                    
                    # Parse header
                    year = int(line[13:17])
                    month = int(line[18:20])
                    day = int(line[21:23])
                    hour = int(line[24:26])
                    reltime = line[27:31].strip()
                    numlev = int(line[32:36])
                    p_src = line[37:45].strip()
                    np_src = line[46:54].strip()
                    
                    # Store date as YYYYMMDD
                    dates[sounding_idx] = year * 10000 + month * 100 + day
                    times[sounding_idx] = hour
                    numlevs[sounding_idx] = numlev
                    
                    # Ensure the entire string is assigned, not just the first character
                    if reltime != '':
                        reltimes.append(reltime)
                    else:
                        reltimes.append('--')
                    if p_src != '':
                        p_srcs.append(p_src)
                    else:
                        p_srcs.append('--')

                    if np_src != '':
                        np_srcs.append(np_src)
                    else:
                        np_srcs.append('--')
                    
                else:
                    # Parse data line
                    try:
                        # IGRA data format: columns are fixed width
                        lvltyp1_val = line[0:1]
                        lvltyp2_val = line[1:2]
                        etime_val = line[3:8].strip()
                        press_val = line[9:15].strip()
                        pflag_val = line[15:16]
                        gph_val = line[16:21].strip()
                        zflag_val = line[21:22]
                        temp_val = line[22:27].strip()
                        tflag_val = line[27:28]
                        rh_val = line[28:33].strip()
                        dpdp_val = line[34:39].strip()
                        wdir_val = line[40:45].strip()
                        wspd_val = line[46:51].strip()
                        
                        # Convert values to float if not empty, otherwise keep as NaN
                        if press_val and press_val != '-9999':
                            pressure[sounding_idx, level_idx] = float(press_val)/100
                        elif press_val:
                            pressure[sounding_idx, level_idx] = float(press_val)
                        
                        if gph_val and gph_val != '-9999' and gph_val != '-8888':
                            gph[sounding_idx, level_idx] = float(gph_val)
                            
                        if temp_val and temp_val != '-9999' and temp_val != '-8888':
                            temp[sounding_idx, level_idx] = float(temp_val)/10
                        elif temp_val:
                            temp[sounding_idx, level_idx] = float(temp_val)
                            
                        if rh_val and rh_val != '-9999' and rh_val != '-8888':
                            rh[sounding_idx, level_idx] = float(rh_val)/10
                        elif rh_val:
                            rh[sounding_idx, level_idx] = float(rh_val)
                            
                        if wdir_val and wdir_val != '-9999' and wdir_val != '-8888':
                            wdir[sounding_idx, level_idx] = float(wdir_val)
                            
                        if wspd_val and wspd_val != '-9999' and wspd_val != '-8888':
                            wspd[sounding_idx, level_idx] = float(wspd_val)/10
                        elif wspd_val:
                            wspd[sounding_idx, level_idx] = float(wspd_val)
                        
                        # Store level type and flags as numeric values if possible
                        if lvltyp1_val and lvltyp1_val.isdigit():
                            lvltyp1[sounding_idx, level_idx] = float(lvltyp1_val)
                            
                        if lvltyp2_val and lvltyp2_val.isdigit():
                            lvltyp2[sounding_idx, level_idx] = float(lvltyp2_val)
                            
                        if etime_val:
                            etime[sounding_idx, level_idx] = float(etime_val)
                            
                        if pflag_val and pflag_val.isalnum():
                            # Store as ASCII value if it's a character
                            pflag[sounding_idx, level_idx] = ord(pflag_val) if len(pflag_val) == 1 else np.nan
                            
                        if zflag_val and zflag_val.isalnum():
                            zflag[sounding_idx, level_idx] = ord(zflag_val) if len(zflag_val) == 1 else np.nan
                            
                        if tflag_val and tflag_val.isalnum():
                            tflag[sounding_idx, level_idx] = ord(tflag_val) if len(tflag_val) == 1 else np.nan
                            
                        if dpdp_val and dpdp_val != '-9999' and dpdp_val != '-8888':
                            dpdp[sounding_idx, level_idx] = float(dpdp_val)/10
                        elif dpdp_val:
                            dpdp[sounding_idx, level_idx] = float(dpdp_val)

                        
                        level_idx += 1
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing line: {line}")
                        print(f"Error details: {e}")
                        # Skip malformed lines
                        continue
    
        
        # Create the NetCDF file
        with nc.Dataset(output_file, 'w', format='NETCDF4') as ncfile:
            # Create dimensions
            ncfile.createDimension('profile', sounding_count)
            ncfile.createDimension('levels', max_levels)
            
            # Create variables
            if main:
                date_var = ncfile.createVariable('date', 'i4', ('profile',))
                time_var = ncfile.createVariable('time', 'i4', ('profile',))
                pressure_var = ncfile.createVariable('pressure', 'f4', ('profile', 'levels'), fill_value=np.nan)
                gph_var = ncfile.createVariable('gph', 'f4', ('profile', 'levels'), fill_value=np.nan)
                temperature_var = ncfile.createVariable('temperature', 'f4', ('profile', 'levels'), fill_value=np.nan)
                rh_var = ncfile.createVariable('rh', 'f4', ('profile', 'levels'), fill_value=np.nan)
                wdir_var = ncfile.createVariable('wdir', 'f4', ('profile', 'levels'), fill_value=np.nan)
                wspd_var = ncfile.createVariable('wspd', 'f4', ('profile', 'levels'), fill_value=np.nan)
                dpdp_var = ncfile.createVariable('dpdp', 'f4', ('profile', 'levels'), fill_value=np.nan)
            else:
                date_var = ncfile.createVariable('date', 'i4', ('profile',))
                time_var = ncfile.createVariable('time', 'i4', ('profile',))
                reltime_var = ncfile.createVariable('reltime', str, ('profile',))
                numlev_var = ncfile.createVariable('numlev', 'i4', ('profile',))
                p_src_var = ncfile.createVariable('p_src', str, ('profile',))
                np_src_var = ncfile.createVariable('np_src', str, ('profile',))
                pressure_var = ncfile.createVariable('pressure', 'f4', ('profile', 'levels'), fill_value=np.nan)
                gph_var = ncfile.createVariable('gph', 'f4', ('profile', 'levels'), fill_value=np.nan)
                temperature_var = ncfile.createVariable('temperature', 'f4', ('profile', 'levels'), fill_value=np.nan)
                rh_var = ncfile.createVariable('rh', 'f4', ('profile', 'levels'), fill_value=np.nan)
                wdir_var = ncfile.createVariable('wdir', 'f4', ('profile', 'levels'), fill_value=np.nan)
                wspd_var = ncfile.createVariable('wspd', 'f4', ('profile', 'levels'), fill_value=np.nan)
                lvltyp1_var = ncfile.createVariable('lvltyp1', 'f4', ('profile', 'levels'), fill_value=np.nan)
                lvltyp2_var = ncfile.createVariable('lvltyp2', 'f4', ('profile', 'levels'), fill_value=np.nan)
                etime_var = ncfile.createVariable('etime', 'f4', ('profile', 'levels'), fill_value=np.nan)
                pflag_var = ncfile.createVariable('pflag', 'f4', ('profile', 'levels'), fill_value=np.nan)
                zflag_var = ncfile.createVariable('zflag', 'f4', ('profile', 'levels'), fill_value=np.nan)
                tflag_var = ncfile.createVariable('tflag', 'f4', ('profile', 'levels'), fill_value=np.nan)
                dpdp_var = ncfile.createVariable('dpdp', 'f4', ('profile', 'levels'), fill_value=np.nan)
            
            # Add attributes
            date_var.units = 'YYYYMMDD'
            date_var.long_name = 'Date of sounding'
            
            time_var.units = 'hour'
            time_var.long_name = 'Time of sounding (UTC)'
            
            pressure_var.units = 'hPa'
            pressure_var.long_name = 'Pressure'
            
            gph_var.units = 'm'
            gph_var.long_name = 'Geopotential height'
            
            # Global attributes
            ncfile.description = f'IGRA sounding data for station {station_name}'
            ncfile.history = f'Created {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
            ncfile.source = f'IGRA v2 data from {data_file}'
            
            # Write data
            if main:
                date_var[:] = dates
                time_var[:] = times
                pressure_var[:, :] = pressure
                gph_var[:, :] = gph
                temperature_var[:, :] = temp
                rh_var[:, :] = rh
                wdir_var[:, :] = wdir
                wspd_var[:, :] = wspd
                dpdp_var[:, :] = dpdp
            else:
                date_var[:] = dates
                time_var[:] = times
                reltime_var[:] = np.array(reltimes)
                numlev_var[:] = numlevs
                p_src_var[:] = np.array(p_srcs)
                np_src_var[:] = np.array(np_srcs)
                pressure_var[:, :] = pressure
                gph_var[:, :] = gph
                temperature_var[:, :] = temp
                rh_var[:, :] = rh
                wdir_var[:, :] = wdir
                wspd_var[:, :] = wspd
                lvltyp1_var[:, :] = lvltyp1
                lvltyp2_var[:, :] = lvltyp2
                etime_var[:, :] = etime
                pflag_var[:, :] = pflag
                zflag_var[:, :] = zflag
                tflag_var[:, :] = tflag
                dpdp_var[:, :] = dpdp
        
        print(f"Successfully created NetCDF file: {output_file}")
        return output_file
    
    except Exception as e:
        print(f"Error creating NetCDF file for station {station_name}: {e}")
        return None
    


@app.route('/')
def index():
    stations = download_locations()
    return render_template('index.html', stations=stations)


@app.route('/station_info/<station_id>')
def station_info(station_id):
    stations = download_locations()
    station = next((s for s in stations if s['id'] == station_id), None)
    return jsonify(station)

@app.route('/availability/<station_id>')
def availability(station_id):
    availability = get_availability(station_id)
    return jsonify(availability)

@app.route('/export/<station_id>')
def export_station(station_id):
    # Create a temporary directory for processing
    temp_dir = os.path.join(os.getcwd(), 'temp', str(datetime.now().timestamp()))
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    try:
        output_file = netcdf_app(station_id, download=True, download_dir=temp_dir)
        if output_file:
            # For AJAX requests, return JSON
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
                # Return a JSON response with file URL
                file_url = f"/download/{os.path.basename(temp_dir)}/{os.path.basename(output_file)}"
                return jsonify({"status": "success", "message": "NetCDF file created successfully", "file_url": file_url})
            else:
                # For direct browser requests, send the file
                response = send_file(output_file, 
                                    as_attachment=True,
                                    download_name=f"{station_id}-full.nc",
                                    mimetype='application/x-netcdf')
                
                # Clean up will happen after the response is sent
                @response.call_on_close
                def cleanup():
                    if os.path.exists(temp_dir):
                        shutil.rmtree(temp_dir)
                
                return response
        else:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            return jsonify({"status": "error", "message": "Failed to create NetCDF file"}), 500
    except Exception as e:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return jsonify({"status": "error", "message": f"Error: {str(e)}"}), 500

@app.route('/download/<temp_folder>/<filename>')
def download_file(temp_folder, filename):
    temp_dir = os.path.join(os.getcwd(), 'temp', temp_folder)
    try:
        response = send_file(os.path.join(temp_dir, filename), 
                            as_attachment=True,
                            download_name=filename,
                            mimetype='application/x-netcdf')
        
        @response.call_on_close
        def cleanup():
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        
        return response
    except Exception as e:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return jsonify({"status": "error", "message": f"Error downloading file: {str(e)}"}), 500
    
@app.route('/export_main/<station_id>')
def export_station_main(station_id):
    # Create a temporary directory for processing
    temp_dir = os.path.join(os.getcwd(), 'temp', str(datetime.now().timestamp()))
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    try:
        output_file = netcdf_app(station_id, download=True, download_dir=temp_dir, main=True)
        if output_file:
            # For AJAX requests, return JSON
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
                # Return a JSON response with file URL
                file_url = f"/download/{os.path.basename(temp_dir)}/{os.path.basename(output_file)}"
                return jsonify({"status": "success", "message": "NetCDF file created successfully", "file_url": file_url})
            else:
                # For direct browser requests, send the file
                response = send_file(output_file, 
                                    as_attachment=True,
                                    download_name=f"{station_id}-main.nc",
                                    mimetype='application/x-netcdf')
                
                # Clean up will happen after the response is sent
                @response.call_on_close
                def cleanup():
                    if os.path.exists(temp_dir):
                        shutil.rmtree(temp_dir)
                
                return response
        else:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
            return jsonify({"status": "error", "message": "Failed to create NetCDF file"}), 500
    except Exception as e:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return jsonify({"status": "error", "message": f"Error: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(debug=True)

