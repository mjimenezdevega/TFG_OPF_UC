using PyPlot, DataFrames

function graphManager(nN, dLine, solGen, solFlows, num_hours)
    # Determinar el tamaño de la malla (aproximadamente cuadrada)
    grid_size = ceil(Int, sqrt(nN))

    # Generar posiciones de los nodos en una cuadrícula
    all_nodes = [(i, j) for i in 1:grid_size for j in 1:grid_size]

    # Asegurar que no exceda la cantidad de nodos disponibles
    nN = min(nN, length(all_nodes))
    nodos = all_nodes[1:nN]

    # Crear diccionario para asignar ID de nodo a coordenadas
    nodo_pos = Dict(i => nodos[i] for i in 1:nN)

    # Crear gráfico para cada hora
    fig, axs = subplots(ceil(Int, sqrt(num_hours)), ceil(Int, sqrt(num_hours)), figsize=(12, 12))

    # Convertir a 1D para iterar fácilmente
    axs = vec(axs)

    # Flag para agregar la leyenda solo una vez
    legend_added = false

    for t in 1:num_hours
        ax = axs[t]

        # Dibujar líneas de conexión entre nodos
        for i in 1:nrow(dLine)
            if haskey(nodo_pos, dLine.tbus[i]) && haskey(nodo_pos, dLine.fbus[i])
                x_start, y_start = nodo_pos[dLine.tbus[i]]
                x_end, y_end = nodo_pos[dLine.fbus[i]]
                ax.plot([x_start, x_end], [y_start, y_end], "k-", linewidth=1)  # Línea negra
            end
        end

        # Dibujar nodos y números de nodos
        x_coords, y_coords = first.(nodos), last.(nodos)
        ax.scatter(x_coords, y_coords, s=50, c="black")  # Nodos en negro
        
        # Etiquetas con el número de nudo
        for i in 1:nN
            x, y = nodo_pos[i]
            ax.text(x + 0.1, y + 0.1, string(i), fontsize=8, color="black")  # Número del nudo al lado del nudo
        end

        # Dibujar flechas de generación (solo si la potencia no es 0)
        for row in eachrow(solGen)
            if row.hour == t
                bus = row.bus
                powerGen = row.powerGen
                if powerGen != 0 && haskey(nodo_pos, bus)
                    x, y = nodo_pos[bus]
                    ax.annotate(
                        text="$(round(powerGen, digits=2))", 
                        xy=(x, y), 
                        xytext=(x, y + 0.2), 
                        arrowprops=Dict(:arrowstyle => "->", :color => "orange", :lw => 2),  # Usamos Dict
                        color="orange", fontsize=10, label="Thermal Power Generation"  # Añadimos etiqueta para leyenda
                    )
                end
            end
        end

        # Dibujar flechas de flujo de potencia (dependiendo de si el flujo es negativo o no)
        for row in eachrow(solFlows)
            if row.hour == t
                fbus = row.fbus
                tbus = row.tbus
                flow = row.flow
                if haskey(nodo_pos, fbus) && haskey(nodo_pos, tbus)
                    x_start, y_start = nodo_pos[fbus]
                    x_end, y_end = nodo_pos[tbus]
                    
                    # Calcular la longitud total de la línea entre los nodos
                    delta_x = x_end - x_start
                    delta_y = y_end - y_start
                    length = sqrt(delta_x^2 + delta_y^2)

                    # Reducir la longitud de las flechas (por ejemplo, un 10% de la longitud total)
                    scale_factor = 0.1  # 10% de la longitud de la línea
                    arrow_length = length * scale_factor

                    # Normalizar las direcciones
                    unit_vector_x = delta_x / length
                    unit_vector_y = delta_y / length

                    # Desplazar las flechas a lo largo de la línea
                    x_start_arrow = x_start + unit_vector_x * arrow_length
                    y_start_arrow = y_start + unit_vector_y * arrow_length
                    x_end_arrow = x_end - unit_vector_x * arrow_length
                    y_end_arrow = y_end - unit_vector_y * arrow_length

                    # Invertir dirección de la flecha si el flujo es negativo
                    if flow < 0
                        # Flecha en sentido contrario (de tbus a fbus)
                        ax.annotate(
                            text="$(abs(flow))", 
                            xy=(x_end_arrow, y_end_arrow), 
                            xytext=(x_start_arrow, y_start_arrow), 
                            arrowprops=Dict(:arrowstyle => "->", :color => "blue", :lw => 2),
                            color="blue", fontsize=8, label="Line Flow"
                        )
                    else
                        # Flecha en dirección normal (de fbus a tbus)
                        ax.annotate(
                            text="$(flow)", 
                            xy=(x_start_arrow, y_start_arrow), 
                            xytext=(x_end_arrow, y_end_arrow), 
                            arrowprops=Dict(:arrowstyle => "->", :color => "blue", :lw => 2),
                            color="blue", fontsize=8, label="Line Flow"
                        )
                    end

                    # Calcular el punto medio de la línea
                    x_mid = (x_start_arrow + x_end_arrow) / 2
                    y_mid = (y_start_arrow + y_end_arrow) / 2

                    # Calcular el ángulo de rotación para alinear el texto con la línea
                    angle = atan(delta_y / delta_x) * (180 / π)

                    # Colocar el texto paralelo a la línea en el punto medio, un poco por encima de la flecha
                    ax.text(
                        x_mid, y_mid + 0.05,  # Un poco arriba de la flecha
                        "$(round(flow, digits=2))", 
                        color="blue", fontsize=6,
                        rotation=angle, ha="center", va="center"
                    )
                end
            end
        end

        # Ajustes de estilo
        ax.set_title("Hora: $t")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(false)

        # Añadir leyenda solo en la primera subgráfica
        if !legend_added
            ax.legend(loc="upper right", fontsize=10)
            legend_added = true
        end
    end

    show()
end