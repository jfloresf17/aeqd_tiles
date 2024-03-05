## Proyecto de Descarga de Datos Geoespaciales

Este proyecto consiste en un script en Python diseñado para descargar datos geoespaciales de diversas fuentes. El script utiliza una función `download` ubicada en el módulo `utils` para descargar los datos especificados.

### Uso

Para utilizar el script, sigue los siguientes pasos:

1. Clona este repositorio en tu máquina local.
2. Asegúrate de tener Python instalado en tu sistema.
3. Instala las dependencias del proyecto. Puedes hacerlo ejecutando `pip install -r requirements.txt`.
4. Ejecuta el script principal `main.py`.

#### Ejemplo de Uso

```python
from utils import download

# Lista de productos de muestra
sample_products = [
    "MODIS_BRDF",  # fecha y bandas
    "MODIS_LAND_SURFACE",  # fecha y bandas
    "GLO_DEM",  # sin fecha
]

# Descargar los productos de muestra para la zona SA y el azulejo T1 especificado
download(sample_products, zone="SA", T1_tile="E0589N0557T1")
```
### Productos Disponibles

El script admite la descarga de los siguientes productos:

- **MODIS_BRDF**: Datos con fecha y bandas.
- **MODIS_LAND_SURFACE**: Datos con fecha y bandas.
- **GLO_DEM**: Datos sin fecha.

(Se pueden agregar más detalles sobre otros productos compatibles).

### Contribución

Si deseas contribuir al proyecto, ¡estamos abiertos a solicitudes de extracción! Si encuentras un error o tienes una idea para una nueva característica, no dudes en abrir un problema.

### Licencia

Este proyecto está bajo la Licencia MIT.

Este README proporciona una breve descripción del proyecto, instrucciones de uso, ejemplos, lista de productos disponibles y detalles sobre cómo contribuir. Recuerda ajustarlo según las necesidades y especificaciones de tu proyecto.
