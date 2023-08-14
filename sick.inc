define sicklaser ranger
(
  sensor
  (
    range [0 8]
    fov 30
    samples 30
  )
  # generic model properties
  size [0.15 0.15 0.15]   # Increase the size for visibility
  color "red"            # Change color to red for better visibility
  show_laser 1           # Makes the laser scan visible
)

