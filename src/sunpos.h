

sky_pos sunpos(time_t t, double lat, double lon); // lat & lon in radians!
void SolarTimes(time_t t, double lat, double lon, time_t * trise, time_t *tnoon, time_t *tset, sky_pos *prise, sky_pos *pnoon, sky_pos *pset);
