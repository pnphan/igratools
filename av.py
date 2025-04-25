import numpy as np
import netCDF4 as nc
import os
import requests
import zipfile
from urllib.parse import urljoin
from datetime import datetime
import io
import shutil
import argparse
import json


def get_availability(station_id, download=True, download_dir=None, download_availability=True):
    """
    Get the availability of IGRA data for a given station.
    """
    availability = []
    try:
                    
        if download_dir is None:
            download_dir = os.path.join(os.getcwd(), str(datetime.now().strftime("%Y-%m-%d")))

        # Create the directory if it doesn't exist
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

        # Create the NetCDF file
        output_file = os.path.join(download_dir, f"{station_id}.nc")


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
                
                # Extract only the text file directly from memory to disk
                print(f"Extracting data for station {station_id}...")
                with zipfile.ZipFile(zip_buffer) as zip_ref:
                    # Get the name of the text file in the zip
                    txt_file_name = None
                    for file_info in zip_ref.infolist():
                        if file_info.filename.endswith('.txt'):
                            txt_file_name = file_info.filename
                            break
                    
                    if txt_file_name:
                        # Extract only the text file
                        with zip_ref.open(txt_file_name) as source, open(os.path.join(download_dir, f"{station_id}-data.txt"), 'wb') as target:
                            shutil.copyfileobj(source, target)
                    else:
                        raise Exception("No text file found in the zip archive")
                
                print(f"Successfully downloaded and extracted data for {station_id} to {download_dir}")
            
            except requests.exceptions.RequestException as e:
                print(f"Error downloading data for station {station_id}: {e}")
                return None
            except zipfile.BadZipFile:
                print(f"Error: The downloaded file for {station_id} is not a valid zip file.")
                return None
            except Exception as e:
                print(f"Unexpected error processing data for {station_id}: {e}")
                return None

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
                else:
                    continue
        availability_grid = np.array(availability)
        years = np.unique(availability_grid[:, 0])

        year_soundings = {}
        for year, count in zip(np.unique(availability_grid[:, 0]), np.bincount(availability_grid[:, 0], minlength=int(np.max(availability_grid[:, 0])+1))[int(np.min(availability_grid[:, 0])):]):
            year_soundings[int(year)] = int(count)

        year_months = {}
        # Count unique months per year
        for year in np.unique(availability_grid[:, 0]):
        # Filter data for this year and get unique months
            year_data = availability_grid[availability_grid[:, 0] == year]
            unique_months = np.unique(year_data[:, 1])
            month_count = len(unique_months)
            year_months[int(year)] = int(month_count)

        year_days = {}
        # Count unique days per year
        for year in np.unique(availability_grid[:, 0]):
            # Filter data for this year and get unique days
            year_data = availability_grid[availability_grid[:, 0] == year]
            # Create a unique identifier for each day of the year (month*100 + day)
            day_of_year = year_data[:, 1] * 100 + year_data[:, 2]
            unique_days = np.unique(day_of_year)
            day_count = len(unique_days)
            year_days[int(year)] = int(day_count)


        if download_availability:
            # Save availability data to a file in the download directory
            availability_file = os.path.join(download_dir, f"{station_id}-availability.json")
            availability_data = {
                'station_id': station_id,
                'raw_data' : availability,
                'num_total_soundings': len(availability),
                'available_years': years.tolist(),
                'num_soundings_per_year': year_soundings,
                'num_months_per_year': year_months,
                'num_days_per_year': year_days
            }
            with open(availability_file, 'w') as f:
                json.dump(availability_data, f, indent=4)
            print(f"Availability data saved to {availability_file}")
            
        return availability_data
    
    except Exception as e:
        print(f"Error fetching availability for station {station_id}: {e}")
        return None
    

def main():
    parser = argparse.ArgumentParser(description='Create NetCDF files from IGRA data for specified stations.')
    parser.add_argument('stations', help='IGRA station identifiers (e.g., USM00072401 USM00072421)') # nargs='+'
    parser.add_argument('--nd', dest='download', action='store_false', 
                        help='Do not download IGRA data files (use existing files)')
    parser.add_argument('--dir', dest='download_dir', type=str, default=None,
                        help='Directory to save IGRA data files (default: current date)')
    parser.add_argument('--da', dest='download_availability', action='store_false', 
                        help='Do not download availability data (use existing files)')
    
    args = parser.parse_args()
    
    # Call the netcdf function with the parsed arguments
    get_availability(args.stations, download=args.download, download_dir=args.download_dir, download_availability=args.download_availability)


if __name__ == "__main__":
    main()