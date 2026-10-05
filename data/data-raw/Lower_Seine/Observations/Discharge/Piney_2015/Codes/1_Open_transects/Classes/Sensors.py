import numpy as np
from Classes.SensorStructure import SensorStructure


class Sensors(object):
    """Class to store data from ADCP sensors.

    Attributes
    ----------
    heading_deg: HeadingData
        Object of HeadingData.
    pitch_deg: SensorStructure
        Pitch data, object of SensorStructure
    roll_deg: SensorStructure
        Roll data, object of SensorStructure
    temperature_deg_c: SensorStructure
        Temperature data, object of SensorStructure
    salinity_ppt: SensorStructure
        Salinity data, object of SensorStructure
    speed_of_sound_mps: SensorStructure
        Speed of sound, object of SensorStructure
    """

    def __init__(self):
        """Initialize class and create variable objects"""

        self.heading_deg = SensorStructure()  # Object of HeadingData
        self.pitch_deg = SensorStructure()  # Pitch data, object of SensorStructure
        self.roll_deg = SensorStructure()  # Roll data, object of SensorStructure
        self.temperature_deg_c = SensorStructure()  # Temperature data, object of SensorStructure
        self.salinity_ppt = SensorStructure()  # Salinity data, object of SensorStructure
        self.speed_of_sound_mps = SensorStructure()  # Speed of sound, object of SensorStructure

    def populate_from_qrev_mat(self, transect):
        """Populates the object using data from previously saved QRev Matlab file.

        Parameters
        ----------
        transect: mat_struct
           Matlab data structure obtained from sio.loadmat
        """
        if hasattr(transect, 'sensors'):
            if hasattr(transect.sensors, 'heading_deg'):
                self.heading_deg.populate_from_qrev_mat(transect.sensors.heading_deg, heading=True)
            if hasattr(transect.sensors, 'pitch_deg'):
                self.pitch_deg.populate_from_qrev_mat(transect.sensors.pitch_deg)
            if hasattr(transect.sensors, 'roll_deg'):
                self.roll_deg.populate_from_qrev_mat(transect.sensors.roll_deg)
            if hasattr(transect.sensors, 'salinity_ppt'):
                self.salinity_ppt.populate_from_qrev_mat(transect.sensors.salinity_ppt)
            if hasattr(transect.sensors, 'speedOfSound_mps'):
                self.speed_of_sound_mps.populate_from_qrev_mat(transect.sensors.speedOfSound_mps)
            if hasattr(transect.sensors, 'temperature_degC'):
                self.temperature_deg_c.populate_from_qrev_mat(transect.sensors.temperature_degC)

    @staticmethod
    def speed_of_sound(temperature, salinity):
        """Computes speed of sound from temperature and salinity.

        Parameters
        ----------
        temperature: float or np.array(float)
            Water temperature at transducer face, in degrees C.
        salinity: float or np.array(float)
            Water salinity at transducer face, in ppt.
        """

        # Not provided in RS Matlab file computed from equation used in TRDI BBSS, from Urick (1983)
        # May not be the same equation as used by SonTek
        sos = 1449.2 + 4.6 * temperature - 0.055 * temperature**2 + 0.00029 * temperature**3 \
            + (1.34 - 0.01 * temperature) * (salinity - 35.0)

        return sos

    @staticmethod
    def unesco_speed_of_sound(t, s, p=10):

        p = p / 10
        sr = np.sqrt(np.abs(s))

        # S ** 2 TERM
        d = 1.727E-3 - 7.9836E-6 * p

        # S ** 3 / 2  TERM
        b1 = 7.3637E-5 + 1.7945E-7 * t
        b0 = -1.922E-2 - 4.42E-5 * t
        b = b0 + b1 * p

        # S ** 1 TERM
        a3 = (-3.389E-13 * t + 6.649E-12) * t + 1.100E-10
        a2 = ((7.988E-12 * t - 1.6002E-10) * t + 9.1041E-9) * t - 3.9064E-7
        a1 = (((-2.0122E-10 * t + 1.0507E-8) * t - 6.4885E-8) * t - 1.2580E-5) * t + 9.4742E-5
        a0 = (((-3.21E-8 * t + 2.006E-6) * t + 7.164E-5) * t - 1.262E-2) * t + 1.389
        a = ((a3 * p + a2) * p + a1) * p + a0

        # S ** 0 TERM
        c3 = (-2.3643E-12 * t + 3.8504E-10) * t - 9.7729E-9
        c2 = (((1.0405E-12 * t - 2.5335E-10) * t + 2.5974E-8) * t - 1.7107E-6) * t + 3.1260E-5
        c1 = (((-6.1185E-10 * t + 1.3621E-7) * t - 8.1788E-6) * t + 6.8982E-4) * t + 0.153563
        c0 = ((((3.1464E-9 * t - 1.47800E-6) * t + 3.3420E-4) * t - 5.80852E-2) * t + 5.03711) * t + 1402.388
        c = ((c3 * p + c2) * p + c1) * p + c0

        # SOUND  SPEED
        sos = c + (a + b * sr + d * s) * s

        return sos

    @staticmethod
    def avg_temperature(transects):
        """Compute mean temperature from temperature data from all transects.

        Parameters
        ----------
        transects: list
            List of TransectData objects
        """

        temps = np.array([])
        for transect in transects:
            if transect.checked:
                temps = np.append(temps, transect.sensors.temperature_deg_c.internal.data)
        return np.nanmean(temps)
