package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

/**
 *  Represents an individual under study.
 */
public final class Sample implements Comparable<Sample> {
    private final String familyID;
    private final String paternalID;
    private final String maternalID;
    private final Sex gender;
    private final String ID;
    private final Affection affection;

    public Sample(final String ID,
                  final String familyID,
                  final String paternalID,
                  final String maternalID,
                  final Sex gender,
                  final Affection affection) {
        Utils.nonNull(ID, "ID is null");
        Utils.nonNull(gender, "sex is null");
        this.ID = ID;
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
        this.affection = affection;
    }

    public Sample(final String ID,
                  final String familyID,
                  final String paternalID,
                  final String maternalID,
                  final Sex gender) {
        this(ID, familyID, paternalID, maternalID, gender, Affection.UNKNOWN);
    }

    public String getID() {
        return ID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public String getPaternalID() {
        return paternalID;
    }

    public String getMaternalID() {
        return maternalID;
    }

    public Sex getSex() {
        return gender;
    }

    public Affection getAffection() {
        return affection;
    }

    @Override
    public String toString() {
        return String.format("Sample %s fam=%s dad=%s mom=%s gender=%s affection=%s", ID, familyID, paternalID, maternalID, gender, affection);
    }

    @Override
    public int compareTo(final Sample sample) {
        return ID.compareTo(sample.getID());
    }

    @Override
    public int hashCode() {
        return ID.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if (o == null) {
            return false;
        }
        if (o instanceof Sample) {
            Sample otherSample = (Sample) o;
            return ID.equals(otherSample.ID) &&
                    equalOrNull(familyID, otherSample.familyID) &&
                    equalOrNull(paternalID, otherSample.paternalID) &&
                    equalOrNull(maternalID, otherSample.maternalID) &&
                    equalOrNull(gender, otherSample.gender) &&
                    equalOrNull(affection, otherSample.affection);
        }
        return false;
    }

    final Sample mergeSamples(final Sample newSample) {
        if (this.equals(newSample))
            return newSample;
        else {
            return new Sample(this.getID(),
                    mergeValues(this.getID(), "FamilyID", this.getFamilyID(), newSample.getFamilyID(), null),
                    mergeValues(this.getID(), "PaternalID", this.getPaternalID(), newSample.getPaternalID(), null),
                    mergeValues(this.getID(), "MaterialID", this.getMaternalID(), newSample.getMaternalID(), null),
                    mergeValues(this.getID(), "Gender", this.getSex(), newSample.getSex(), Sex.UNKNOWN),
                    mergeValues(this.getID(), "Affection", this.getAffection(), newSample.getAffection(), Affection.UNKNOWN)
            );
        }
    }

    private static <T> T mergeValues(final String name, final String field, final T o1, final T o2, final T emptyValue) {
        if (o1 == null || o1.equals(emptyValue)) {
            // take o2 if both are null, otherwise keep o2
            return o2 == null ? null : o2;
        }
        else {
            if (o2 == null || o2.equals(emptyValue)) {
                return o1; // keep o1, since it's a real value
            }
            else { // both o1 and o2 have a value
                if (o1 instanceof String && o1.equals(o2)) {
                    return o1;
                }
                else if (o1 == o2) {
                    return o1;
                }
                else {
                    throw new UserException("Inconsistent values detected for " + name + " for field " + field + " value1 " + o1 + " value2 " + o2);
                }
            }
        }
    }

    private static boolean equalOrNull(final Object o1, final Object o2) {
        if (o1 == null) {
            return o2 == null;
        }
        else {
            return o2 == null ? false : o1.equals(o2);
        }
    }

}
