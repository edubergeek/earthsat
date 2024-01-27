from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
from astropy.coordinates import ITRS, ICRS, GCRS
from astropy import units as u
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS
import astropy
import sgp4

from platform import python_version

class TLE:
  def PrintVersion(self):
    print("AstroPy: ", astropy.__version__)
    print("SGP4:    ", sgp4.__version__)
    print("Python:  ", python_version())

  def GetSites():
    return EarthLocation.get_site_names()
  
  def TEME(self, tle, obstime):
    if len(tle) == 2:
      self.sat = Satrec.twoline2rv(tle[0], tle[1])
    else:
      return None

    # Use satellite.sgp4 to compute the TEME position and velocity at observation time:
    error_code, teme_p, teme_v = self.sat.sgp4(obstime.jd1, obstime.jd2)  # in km and km/s
    if error_code != 0:
      raise RuntimeError(SGP4_ERRORS[error_code])

    # Cache observing time
    self.obstime = obstime

    # Get the position in km and veloctiy in km/sec
    self.teme_pos = CartesianRepresentation(teme_p*u.km)
    self.teme_vel = CartesianDifferential(teme_v*u.km/u.s)

    # Create a position in the TEME reference frame at the time of observation:
    self.TEME = TEME(self.teme_pos.with_differentials(self.teme_vel), obstime=obstime)
    return self.TEME
  
  def SatPos(self, obslocation, teme=None, obstime=None):
    # Check for use of cached TEME and obstime
    if teme is None:
      teme = self.TEME
    if obstime is None:
      obstime = self.obstime

    # Find the overhead latitude, longitude, and height of the satellite:
    itrs = ITRS(obstime=obstime)
    itrs_geo = teme.transform_to(itrs)
    satlocation = itrs_geo.earth_location

    # Cache the location
    self.obslocation = obslocation
    self.satlocation = satlocation
  
    # Translate to the latitude, longitude, and height of the satellite as seen from the location:
    topo_itrs_repr = itrs_geo.cartesian.without_differentials() - obslocation.get_itrs(obstime).cartesian
    itrs = ITRS(topo_itrs_repr, obstime=obstime, location=obslocation)
    self.ITRS = itrs

    # Transform to terrestrial coords (altitude, azimuth)
    aa = itrs.transform_to(AltAz(obstime=obstime, location=obslocation))
  
    # Transform to celestial coords (right ascension, declination)
    icrs = itrs.transform_to(ICRS())
    self.ICRS = icrs
  
    gcrs = obslocation.get_gcrs(obstime)
    self.GCRS = gcrs

    self.altaz = AltAz(alt=aa.alt, az=aa.az, location=obslocation, obstime=obstime)
    
    # Convert to equatorial coordinates
    self.radec = self.altaz.transform_to(gcrs)
    #print(eqc.ra*'deg', eqc.dec*'deg')
  
    itrs_geo_p = itrs_geo.cartesian.without_differentials()
    itrs_geo_v = itrs_geo.cartesian.differentials['s']
  
    topo_itrs_p = itrs_geo_p - obslocation.get_itrs(obstime).cartesian
    topo_itrs_repr = topo_itrs_p.with_differentials(itrs_geo_v)
  
    itrs_v = ITRS(topo_itrs_repr, obstime=obstime, location=obslocation)
    #print(itrs_v)
    icrs_v = itrs_v.transform_to(ICRS())
    #print(icrs_v)
    altaz_v = itrs_v.transform_to(self.altaz)

    self.alt = self.altaz.alt * 'deg'
    self.az = self.altaz.az * 'deg'
    self.vx = itrs_v.v_x * 'km/s'
    self.vy = itrs_v.v_y * 'km/s'
    self.vz = itrs_v.v_z * 'km/s'
    self.ra = self.radec.ra * 'deg'
    self.dec = self.radec.dec * 'deg'
    self.pm_ra = icrs_v.pm_ra_cosdec *1.0*u.year / (365.25*86400*u.second)
    self.pm_dec = icrs_v.pm_dec *1.0*u.year / (365.25*86400*u.second)
    self.pm_az = altaz_v.pm_az_cosalt *1.0*u.year / (365.25*86400*u.second)
    self.pm_alt = altaz_v.pm_alt *1.0*u.year / (365.25*86400*u.second)
    self.rv = icrs_v.radial_velocity * 'km/s'



