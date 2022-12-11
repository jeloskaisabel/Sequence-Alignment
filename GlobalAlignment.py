# %% Imports
from numpy import full
import pandas as pd
from Alignment import Alignment

class GlobalAlignment(object):
    def __init__(self, seq1, seq2, match_score, mismatch_penalty, gap_penalty):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match_bonus = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
        n_rows = len("-" + self.seq1)
        n_columns = len("-" + self.seq2)
        self.scoring_array = full([n_rows, n_columns], 0.0)
        self.traceback_array = full([n_rows, n_columns], "-")
        self.row_labels = [label for label in "-" + self.seq1]
        self.column_labels = [label for label in "-" + self.seq2]
        self.ans = []
        self.dfScoringArray = pd.DataFrame(self.scoring_array, index=self.row_labels, columns=self.column_labels)
        self.dfTracebackArray = pd.DataFrame(self.traceback_array, index=self.row_labels, columns=self.column_labels)
        print("Scoring array:\n", self.scoring_array)
        print("Traceback array:\n", self.traceback_array)

    def globalAlignment(self):
        # %% Implementando el algoritmo de Needleman-Wunsch
        # construir una matriz de ceros
        n_rows = len(self.seq1) + 1  # una fila extra arriba
        n_columns = len(self.seq2) + 1  # una columna extra a la izquierda

        self.scoring_array = full([n_rows, n_columns], 0.0)
        self.traceback_array = full([n_rows, n_columns], "-")

        # Definimos las flechas Unicode que usaremos en la matriz de seguimiento
        up_arrow = "\u2191"
        left_arrow = "\u2190"
        up_left_arrow = "\u2196"

        arrow = "-"

        # iterar sobre las columnas primero porque queremos hacer
        # todas las columnas para la fila 1 antes de la fila 2
        for row in range(n_rows):
            for col in range(n_columns):
                if row == 0 and col == 0:
                    # Estamos en la esquina superior derecha
                    score = 0
                    arrow = "-"
                elif row == 0:
                    # Estamos en la primera fila pero NO en la esquina

                    # Buscamos la puntuación de la celda anterior
                    # (a la izquierda) en la matriz de puntuación\
                    previous_score = self.scoring_array[row, col - 1]
                    # añadimos la penalización por brecha a la puntuación
                    score = previous_score + self.gap_penalty
                    arrow = left_arrow
                elif col == 0:
                    # Estamos en la primera columna pero no en la primera fila
                    previous_score = self.scoring_array[row - 1, col]
                    score = previous_score + self.gap_penalty
                    arrow = up_arrow
                else:
                    # Estamos en una celda 'media' de la alineación

                    # Calcule los puntajes por venir desde arriba,
                    # desde la izquierda (lo que representa una inserción en seq1)
                    cell_to_the_left = self.scoring_array[row, col - 1]
                    from_left_score = cell_to_the_left + self.gap_penalty

                    # o desde arriba (que representa una inserción en seq2)
                    above_cell = self.scoring_array[row - 1, col]
                    from_above_score = above_cell + self.gap_penalty

                    # celda diagonal, que representa una sustitución (por ejemplo, A --> T)
                    diagonal_left_cell = self.scoring_array[row - 1, col - 1]

                    # NOTA: dado que la tabla tiene una fila y una columna adicionales
                    # (las que están en blanco), al indexar de nuevo a la secuencia queremos
                    # fila -1 y columna - 1. ya que la fila 1 representa el carácter 0 de la secuencia.
                    if self.seq1[row - 1] == self.seq2[col - 1]:
                        diagonal_left_cell_score = diagonal_left_cell + self.match_bonus
                    else:
                        diagonal_left_cell_score = diagonal_left_cell + self.mismatch_penalty

                    score = max([from_left_score, from_above_score, diagonal_left_cell_score])
                    # tomamos el máximo

                    # notese de qué celda era la máxima en la matriz de rastreo
                    # usando flechas Unicode
                    if score == from_left_score:
                        arrow = left_arrow
                    elif score == from_above_score:
                        arrow = up_arrow
                    elif score == diagonal_left_cell_score:
                        arrow = up_left_arrow

                self.traceback_array[row, col] = arrow
                self.scoring_array[row, col] = score

        # %% Traceback
        '''
        Alineamos seq1 y seq2 usando la matriz de rastreo y regrese como dos cadenas
        - traceback_array: una matriz numpy con caracteres de flecha que indican la dirección 
        desde la que se originó el mejor camino a una posición de alineación determinada
        - seq1 - una secuencia representada como una cadena
        - seq2 - una secuencia representada como una cadena
        - up_arrow - el Unicode utilizado para las flechas hacia arriba (hay varios símbolos de flecha en Unicode)
        - left_arrow - el Unicode usado para las flechas izquierdas
        - up_left_arrow - el Unicode utilizado para las flechas diagonales
        - stop - el símbolo utilizado en la parte superior izquierda para indicar el final de la alineación
        '''
        up_arrow = "\u2191"
        left_arrow = "\u2190"
        up_left_arrow = "\u2196"
        stop = "-"

        row = len(self.seq1)
        col = len(self.seq2)
        arrow = self.traceback_array[row, col]
        aligned_seq1 = ""
        aligned_seq2 = ""
        alignment_indicator = ""

        while arrow != "-":
            self.ans.append("Fila actual: " + str(row)+"\n")
            print("Fila actual:", row)
            self.ans.append("Columna actual: " + str(col) + "\n")
            print("Columna actual:", col)
            arrow = self.traceback_array[row, col]
            self.ans.append("Flecha: " + arrow + "\n")
            print("Flecha:", arrow)

            if arrow == up_arrow:
                self.ans.append("Insertar indel en la secuencia superior\n")
                print("Insertar indel en la secuencia superior")
                # We want to add the new indel onto the left
                # side of the growing aligned sequence
                aligned_seq2 = "-" + aligned_seq2
                aligned_seq1 = self.seq1[row - 1] + aligned_seq1
                alignment_indicator = " " + alignment_indicator
                row -= 1

            elif arrow == up_left_arrow:
                self.ans.append("Coincidencia o discrepancia\n")
                print("Coincidencia o discrepancia")
                # Note that we look up the row-1 and col-1 indexes
                # because there is an extra "-" character at the
                # start of each sequence
                seq1_character = self.seq1[row - 1]
                seq2_character = self.seq2[col - 1]
                aligned_seq1 = self.seq1[row - 1] + aligned_seq1
                aligned_seq2 = self.seq2[col - 1] + aligned_seq2
                if seq1_character == seq2_character:
                    alignment_indicator = "|" + alignment_indicator
                else:
                    alignment_indicator = " " + alignment_indicator
                row -= 1
                col -= 1

            elif arrow == left_arrow:
                self.ans.append("Insertar indel en la secuencia izquierda\n")
                print("Insertar indel en la secuencia izquierda")
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = self.seq2[col - 1] + aligned_seq2
                alignment_indicator = " " + alignment_indicator
                col -= 1

            elif arrow == stop:
                break
            else:
                raise ValueError(
                    f"Traceback array entry at {row},{col}: {arrow} is not recognized as an up arrow ({up_arrow}),left_arrow ({left_arrow}), up_left_arrow ({up_left_arrow}), or a stop ({stop}).")

        self.dfScoringArray = pd.DataFrame(self.scoring_array, index=self.row_labels, columns=self.column_labels)
        self.dfTracebackArray = pd.DataFrame(self.traceback_array, index=self.row_labels, columns=self.column_labels)
        ans = Alignment(aligned_seq1, aligned_seq2, self.scoring_array[len(self.seq1)][len(self.seq2)])
        return ans

