#!/usr/bin/env python3
"""
NetCDF MCP Server for Claude Code
Provides tools for reading and analyzing NetCDF files using the FastMCP framework.
"""

import os
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any
from fastmcp import FastMCP

# Initialize the MCP server
mcp = FastMCP("NetCDF Tools")

# Set the base directory for NetCDF files
NETCDF_BASE_DIR = os.environ.get("NETCDF_BASE_DIR", ".")

def _validate_file_path(file_path: str) -> str:
    """Validate and resolve file path within the base directory."""
    base_path = Path(NETCDF_BASE_DIR).resolve()
    full_path = (base_path / file_path).resolve()
    
    # Ensure the file is within the base directory (security check)
    if not str(full_path).startswith(str(base_path)):
        raise ValueError(f"File path must be within {base_path}")
    
    if not full_path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {full_path}")
    
    return str(full_path)

@mcp.tool
def list_netcdf_files(directory: str = ".") -> List[str]:
    """
    List all NetCDF files in the specified directory.
    
    Args:
        directory: Directory to search (relative to base directory)
    
    Returns:
        List of NetCDF file paths
    """
    try:
        base_path = Path(NETCDF_BASE_DIR).resolve()
        search_path = (base_path / directory).resolve()
        
        # Security check
        if not str(search_path).startswith(str(base_path)):
            raise ValueError(f"Directory must be within {base_path}")
        
        netcdf_files = []
        for ext in ["*.nc", "*.nc4", "*.netcdf"]:
            netcdf_files.extend(search_path.glob(ext))
        
        # Return relative paths from base directory
        return [str(f.relative_to(base_path)) for f in netcdf_files]
    
    except Exception as e:
        return [f"Error listing files: {str(e)}"]

@mcp.tool
def ncdump_header(file_path: str) -> str:
    """
    Get NetCDF file header information (dimensions, variables, attributes).
    Equivalent to 'ncdump -h filename.nc'
    
    Args:
        file_path: Path to the NetCDF file (relative to base directory)
    
    Returns:
        Header information as text
    """
    try:
        full_path = _validate_file_path(file_path)
        result = subprocess.run(
            ["ncdump", "-h", full_path],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            return f"Error running ncdump: {result.stderr}"
        
        return result.stdout
    
    except subprocess.TimeoutExpired:
        return "Error: ncdump command timed out"
    except FileNotFoundError:
        return "Error: ncdump command not found. Please ensure NetCDF tools are installed."
    except Exception as e:
        return f"Error: {str(e)}"

@mcp.tool
def ncdump_coords(file_path: str) -> str:
    """
    Get NetCDF coordinate variable data.
    Equivalent to 'ncdump -c filename.nc'
    
    Args:
        file_path: Path to the NetCDF file (relative to base directory)
    
    Returns:
        Coordinate data as text
    """
    try:
        full_path = _validate_file_path(file_path)
        result = subprocess.run(
            ["ncdump", "-c", full_path],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode != 0:
            return f"Error running ncdump: {result.stderr}"
        
        return result.stdout
    
    except subprocess.TimeoutExpired:
        return "Error: ncdump command timed out"
    except Exception as e:
        return f"Error: {str(e)}"

@mcp.tool
def ncdump_variable(file_path: str, variable_name: str) -> str:
    """
    Get data for a specific variable from NetCDF file.
    Equivalent to 'ncdump -v variable_name filename.nc'
    
    Args:
        file_path: Path to the NetCDF file (relative to base directory)
        variable_name: Name of the variable to extract
    
    Returns:
        Variable data as text
    """
    try:
        full_path = _validate_file_path(file_path)
        result = subprocess.run(
            ["ncdump", "-v", variable_name, full_path],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            return f"Error running ncdump: {result.stderr}"
        
        return result.stdout
    
    except subprocess.TimeoutExpired:
        return "Error: ncdump command timed out"
    except Exception as e:
        return f"Error: {str(e)}"

@mcp.tool
def get_netcdf_info(file_path: str) -> Dict[str, Any]:
    """
    Get comprehensive information about a NetCDF file using Python netCDF4 library.
    
    Args:
        file_path: Path to the NetCDF file (relative to base directory)
    
    Returns:
        Dictionary with file information
    """
    try:
        import netCDF4 as nc
        
        full_path = _validate_file_path(file_path)
        
        with nc.Dataset(full_path, 'r') as dataset:
            info = {
                "file_format": dataset.data_model,
                "dimensions": {name: len(dim) for name, dim in dataset.dimensions.items()},
                "variables": {},
                "global_attributes": {attr: getattr(dataset, attr) for attr in dataset.ncattrs()}
            }
            
            # Get variable information
            for var_name, var in dataset.variables.items():
                info["variables"][var_name] = {
                    "dimensions": var.dimensions,
                    "shape": var.shape,
                    "dtype": str(var.dtype),
                    "attributes": {attr: getattr(var, attr) for attr in var.ncattrs()}
                }
        
        return info
    
    except ImportError:
        return {"error": "netCDF4 library not available. Please install: pip install netCDF4"}
    except Exception as e:
        return {"error": str(e)}

@mcp.tool
def read_netcdf_variable_summary(file_path: str, variable_name: str, max_values: int = 10) -> Dict[str, Any]:
    """
    Read a NetCDF variable and provide summary statistics.
    
    Args:
        file_path: Path to the NetCDF file (relative to base directory)
        variable_name: Name of the variable to read
        max_values: Maximum number of values to show in preview
    
    Returns:
        Dictionary with variable summary and sample data
    """
    try:
        import netCDF4 as nc
        import numpy as np
        
        full_path = _validate_file_path(file_path)
        
        with nc.Dataset(full_path, 'r') as dataset:
            if variable_name not in dataset.variables:
                return {"error": f"Variable '{variable_name}' not found in file"}
            
            var = dataset.variables[variable_name]
            data = var[:]
            
            summary = {
                "variable_name": variable_name,
                "shape": var.shape,
                "dimensions": var.dimensions,
                "dtype": str(var.dtype),
                "attributes": {attr: getattr(var, attr) for attr in var.ncattrs()}
            }
            
            # Add statistics for numeric data
            if np.issubdtype(data.dtype, np.number):
                # Handle masked arrays
                if hasattr(data, 'mask'):
                    valid_data = data.compressed()
                else:
                    valid_data = data.flatten()
                
                if len(valid_data) > 0:
                    summary["statistics"] = {
                        "min": float(np.min(valid_data)),
                        "max": float(np.max(valid_data)),
                        "mean": float(np.mean(valid_data)),
                        "std": float(np.std(valid_data)),
                        "valid_points": len(valid_data),
                        "total_points": data.size
                    }
            
            # Add sample data
            flat_data = data.flatten()
            sample_size = min(max_values, len(flat_data))
            summary["sample_data"] = [str(x) for x in flat_data[:sample_size]]
            
            return summary
    
    except ImportError:
        return {"error": "netCDF4 and numpy libraries required. Please install: pip install netCDF4 numpy"}
    except Exception as e:
        return {"error": str(e)}

@mcp.resource("netcdf://help")
def netcdf_help() -> str:
    """Get help information about available NetCDF tools."""
    return """
NetCDF MCP Server - Available Tools:

1. list_netcdf_files(directory="."): List NetCDF files in directory
2. ncdump_header(file_path): Get file structure (ncdump -h)
3. ncdump_coords(file_path): Get coordinate data (ncdump -c)
4. ncdump_variable(file_path, variable_name): Get specific variable data
5. get_netcdf_info(file_path): Get comprehensive file info using Python
6. read_netcdf_variable_summary(file_path, variable_name): Get variable statistics

Environment Variables:
- NETCDF_BASE_DIR: Base directory for NetCDF files (default: current directory)

Usage Examples:
- list_netcdf_files("data/")
- ncdump_header("myfile.nc")
- get_netcdf_info("climate_data.nc4")
- read_netcdf_variable_summary("data.nc", "temperature")
"""

if __name__ == "__main__":
    # Run the MCP server
    mcp.run()